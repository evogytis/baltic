import logging
import warnings
import os
from collections import Counter
import numpy as np
import matplotlib as mpl
from scipy.interpolate import interp1d
from scipy.optimize import minimize
from baltic.bt_utils import calendar_to_decimal_date, desaturate_cmap
from baltic import Reticulation

logger = logging.getLogger("baltic.curonia")

def indexObs(node, count = 0):
    """
    Traverse tree structure, assign an index to every observation.
    """
    # node.traits['tps']=[bt.decimalDate(meta[strain]['date']) for strain in node['strains']]
    descendants = len(node.traits['tps'])
    
    node.traits['leafs'] = np.array(range(count, count + descendants), dtype = np.int_)
    
    count += descendants
    
    if node.is_node():
        for child in node.children:
            count = indexObs(child, count)
    
    return count


def propagateParentObs(node, traitName):
    """
    Propagate observations to parent node from unique counts of each branch.
    """
    if node.is_node():
        for child in node.children:
            node.traits[traitName] = np.hstack((node.traits[traitName], propagateParentObs(child, traitName)))
    
    return node.traits[traitName]



def make_pivots(pivots, tps):
    '''
    if pivots is a scalar, make a grid of pivot points covering the entire range
    Parameters
    ----------
    pivots : scalar or iterable
        either number of pivots (a scalar) or the actual pivots
        (will be cast to array and returned)
    tps : np.array
        observation time points. Will generate pivots spanning min/max
    Returns
    -------
    pivots : np.array
        array of pivot values
    '''
    if np.isscalar(pivots):
        tps = np.asarray(tps)
        if tps.size < 2:
            raise ValueError("tps must contain at least two time points to construct pivots.")  # CHANGED: guard
        dt = np.max(tps) - np.min(tps)
        return np.linspace(np.min(tps) - 0.01 * dt, np.max(tps) + 0.01 * dt, int(pivots))
    else:
        pivots = np.asarray(pivots)
        if pivots.size < 2:
            raise ValueError("Need at least two pivot points.")  # CHANGED: guard
        return pivots


def count_observations(pivots, tps):
    """
    Count how many individual observations fall into each pivot interval.
    """
    pivots = make_pivots(pivots, tps)
    dt = pivots[1] - pivots[0]
    # bins are centered on pivots, each bin spanning +/- 0.5*dt around each pivot
    counts, _ = np.histogram(tps, bins=[pivots[0] - 0.5 * dt] + list(pivots + 0.5 * dt))
    return counts


def running_average(obs, ws):
    '''
    calculates a running average
    obs     --  observations
    ws      --  window size (number of points to average)
    Parameters
    ----------
    obs : list/np.array(bool)
        observations
    ws : int
        window size as measured in number of consecutive points
    Returns
    -------
    np.array(float)
        running average of the boolean observations
    '''
    ws = int(ws)
    try:
        tmp_vals = np.convolve(np.ones(ws, dtype=float) / ws, obs, mode='same')
        # fix the edges. using mode='same' assumes zeros outside the range
        if ws % 2 == 0:
            tmp_vals[:ws // 2] *= float(ws) / np.arange(ws // 2, ws)
            if ws // 2 > 1:
                tmp_vals[-ws // 2 + 1:] *= float(ws) / np.arange(ws - 1, ws // 2, -1.0)
        else:
            tmp_vals[:ws // 2] *= float(ws) / np.arange(ws // 2 + 1, ws)
            tmp_vals[-ws // 2:] *= float(ws) / np.arange(ws, ws // 2, -1.0)
    except Exception:
        tmp_vals = 0.5 * np.ones_like(obs, dtype=float)
    return tmp_vals


def fix_freq(freq, pc):
    '''
    restricts frequencies to the interval [pc, 1-pc]
    removes np.nan values and avoids taking logarithms of 0 or divisions by 0
    '''
    freq = np.asarray(freq)
    freq = freq.copy()
    freq[np.isnan(freq)] = pc
    return np.minimum(1 - pc, np.maximum(pc, freq))


def logit_transform(freq, pc):
    """
    Transform frequencies to logit space after pseudocount clipping.

    Parameters
    ----------
    freq : array-like
        Frequencies to transform.

    pc : float
        Pseudocount used to clip frequencies away from 0 and 1.

    Returns
    -------
    np.ndarray
        Logit-transformed frequencies.
    """
    freq = np.asarray(freq)
    f = fix_freq(freq, pc)
    return np.log(f / (1 - f))


def logit_inv(logit_freq, pc):
    """
    Transform logit-space values back to clipped frequencies.

    Parameters
    ----------
    logit_freq : array-like
        Values in logit space.

    pc : float
        Pseudocount used to clip the reconstructed frequencies.

    Returns
    -------
    np.ndarray
        Frequencies constrained to the interval ``[pc, 1 - pc]``.
    """
    logit_freq = np.asarray(logit_freq)
    tmp = np.exp(logit_freq)
    tmp = np.maximum(pc / (1 - pc), np.minimum((1 - pc) / pc, tmp))  # CHANGED: tighter symmetric clipping
    return tmp / (1.0 + tmp)


def pq(p):
    """
    Compute the Bernoulli variance term ``p * (1 - p)``.

    Parameters
    ----------
    p : array-like
        Frequencies or probabilities.

    Returns
    -------
    np.ndarray
        Elementwise Bernoulli variance term.
    """
    p = np.asarray(p)
    return p * (1 - p)


# ------------------------------------------------------------
# SINGLE-CLADE FREQUENCY ESTIMATOR
# ------------------------------------------------------------

class frequency_estimator(object):
    '''
    estimates a smooth frequency trajectory given a series of time stamped
    0/1 observations. The most likely set of frequencies at specified pivot values
    is determined by numerical minimization. Likelihood consists of a Bernoulli
    sampling term as well as a term penalizing rapid frequency shifts.
    '''

    def __init__(self, tps, obs, pivots,
                 stiffness=20.0, inertia=0.0,
                 tol=1e-3, pc=1e-4, ws=100, log_thres=10,
                 method='powell', **kwargs):
        """
        Parameters
        ----------
        tps : list/np.array(float)
            array with numerical dates (one per sample)
        obs : list/np.array(bool)
            array with boolean observations (same length as tps)
        pivots : int/np.array(float)
            either integer specifying the number of pivot values,
            or list of explicit pivots
        """
        tmp_obs = np.array(sorted(zip(tps, obs), key=lambda x: x[0]))
        self.tps = tmp_obs[:, 0]
        self.obs = np.array(tmp_obs[:, 1], dtype=bool)
        self.stiffness = stiffness

        self.inertia = inertia
        self.interpolation_type = 'linear'
        self.tol = tol
        self.reg = 1e-6
        self.pc = pc
        self.ws = ws
        self.log_thres = log_thres
        self.verbose = 0
        self.method = method

        self.pivots = make_pivots(pivots, self.tps)

        good_tps = (self.tps >= self.pivots[0]) & (self.tps < self.pivots[-1])
        self.tps = self.tps[good_tps]
        self.obs = self.obs[good_tps]

    def initial_guess(self, pc=0.01):
        """
        Construct an initial frequency trajectory estimate on the pivot grid.

        Parameters
        ----------
        pc : float, optional
            Pseudocount used to clip the initial estimate away from 0 and 1.

        Returns
        -------
        np.ndarray
            Initial frequency values at each pivot.
        """
        # generate a useful initial guess from a running average of the counts
        if self.tps.size == 0:
            # CHANGED: guard against empty tps
            return np.full_like(self.pivots, 0.5, dtype=float)

        if self.ws < len(self.obs):
            tmp_vals = running_average(self.obs, self.ws)
        else:
            tmp_vals = running_average(self.obs, len(self.obs))

        tmp_interpolator = interp1d(self.tps, tmp_vals, bounds_error=False, fill_value=-1)
        pivot_freq = tmp_interpolator(self.pivots)
        pivot_freq[self.pivots <= tmp_interpolator.x[0]] = tmp_vals[0]
        pivot_freq[self.pivots >= tmp_interpolator.x[-1]] = tmp_vals[-1]
        pivot_freq = fix_freq(pivot_freq, pc)
        return pivot_freq

    def stiffLH(self):
        """
        Compute the smoothness penalty term of the trajectory likelihood.

        Returns
        -------
        float
            Log-likelihood contribution from the diffusion-style smoothness
            prior.
        """
        freq = self.pivot_freq
        dfreq = np.diff(freq)
        dfreq_penalty = dfreq[1:] - self.inertia * dfreq[:-1]
        pq_dt_freq = self.dt * pq(freq[:-1])

        # return wright-fisher diffusion likelihood for frequency change.
        return -0.25 * self.stiffness * (
            np.sum(dfreq_penalty ** 2 / pq_dt_freq[1:]) + dfreq[0] ** 2 / pq_dt_freq[0]
        )

    def learn(self, initial_guess=None):
        """
        Fit the frequency trajectory by numerical optimization.

        Parameters
        ----------
        initial_guess : callable, optional
            Function that accepts the pivot grid and returns an initial
            frequency trajectory.
        """
        self.dt = np.diff(self.pivots)

        def logLH(x):
            # x is logit frequency
            self.pivot_freq = logit_inv(x, self.pc)
            try:
                freq = interp1d(self.pivots, x, kind=self.interpolation_type,
                                bounds_error=False, assume_sorted=True)
            except Exception:
                freq = interp1d(self.pivots, x, kind=self.interpolation_type,
                                bounds_error=False)
            estfreq = freq(self.tps)
            # Bernoulli log-likelihood over samples: obs*x - log(1+e^x)
            bernoulli_LH = np.sum(estfreq[self.obs]) - np.sum(np.log(1 + np.exp(estfreq)))

            stiffness_LH = self.stiffLH()
            LH = stiffness_LH + bernoulli_LH
            # penalize extreme logits beyond log_thres
            return -LH + np.sum((np.abs(x) > self.log_thres) * x ** 2)

        if initial_guess is None:
            initial_freq = self.initial_guess(pc=0.01)
        else:
            initial_freq = initial_guess(self.pivots)

        self.frequency_estimate = interp1d(self.pivots, initial_freq,
                                           kind=self.interpolation_type,
                                           bounds_error=False)
        self.pivot_freq = self.frequency_estimate(self.pivots)

        self.sol = minimize(logLH, logit_transform(self.pivot_freq, pc=self.pc),
                            method=self.method)
        if self.sol['success']:
            self.pivot_freq = logit_inv(self.sol['x'], self.pc)
        else:
            if self.method != "powell":
                print("Optimization failed, trying with powell")
                self.sol = minimize(logLH, logit_transform(self.pivot_freq, self.pc),
                                    method='powell')
            self.pivot_freq = logit_inv(self.sol['x'], self.pc)

        self.pivot_freq = fix_freq(self.pivot_freq, self.pc)
        self.frequency_estimate = interp1d(self.pivots, self.pivot_freq,
                                           kind=self.interpolation_type,
                                           bounds_error=False)

        if self.verbose:
            print("neg logLH using", len(self.pivots), "pivots:", self.sol['fun'])


class freq_est_clipped(object):
    """
    Wrapper for frequency_estimator that restricts estimation to a sensible
    time window around observations to avoid optimizing where there is no data.
    """

    def __init__(self, tps, obs, pivots, name=None, dtps=None, **kwargs):
        """
        Initialize a clipped frequency estimator around the observation window.

        Parameters
        ----------
        tps : array-like
            Observation times.

        obs : array-like
            Boolean observations aligned with *tps*.

        pivots : np.ndarray
            Pivot grid on which frequencies will be estimated.

        name : str, optional
            Label used in diagnostic messages.

        dtps : float, optional
            Additional margin around the first and last positive observation.

        **kwargs
            Additional keyword arguments forwarded to
            :class:`frequency_estimator`.
        """
        super(freq_est_clipped, self).__init__()
        tmp_obs = np.array(sorted(zip(tps, obs), key=lambda x: x[0]))
        self.tps = tmp_obs[:, 0]
        self.obs = np.array(tmp_obs[:, 1], dtype=bool)
        self.pivots = pivots
        self.name = name

        # pivot spacing
        pivot_dt = self.pivots[1] - self.pivots[0]
        if dtps is None:
            self.dtps = 2.0 * pivot_dt
        else:
            self.dtps = max(dtps, pivot_dt)

        # CHANGED: simple, robust window based on first/last positive obs
        pos_idx = np.where(self.obs)[0]
        if pos_idx.size == 0:
            print(f"[{self.name}] no positive observations")
            self.valid = False
            return

        first_t = self.tps[pos_idx[0]]
        last_t = self.tps[pos_idx[-1]]

        tps_lower_cutoff = max(self.pivots[0], first_t - self.dtps)
        tps_upper_cutoff = min(self.pivots[-1], last_t + self.dtps)

        self.good_tps = (self.tps >= tps_lower_cutoff) & (self.tps <= tps_upper_cutoff)
        self.valid = True
        if self.good_tps.sum() < 3:
            print(f"[{self.name}] too few valid time points:", self.good_tps.sum())
            self.valid = False
            return

        reduced_obs = self.obs[self.good_tps]
        reduced_tps = self.tps[self.good_tps]

        self.pivot_lower_cutoff = min(reduced_tps[0], tps_lower_cutoff) - pivot_dt
        self.pivot_upper_cutoff = max(reduced_tps[-1], tps_upper_cutoff) + pivot_dt

        self.good_pivots = (self.pivots >= self.pivot_lower_cutoff) & \
                           (self.pivots <= self.pivot_upper_cutoff)

        if self.good_pivots.sum() < 2:
            # expand window slightly if needed
            self.good_pivots = self.good_pivots | np.roll(self.good_pivots, 1) | np.roll(self.good_pivots, -1)

        self.fe = frequency_estimator(reduced_tps, reduced_obs,
                                      self.pivots[self.good_pivots], **kwargs)

    def learn(self):
        """
        Fit the clipped frequency estimator and expand the result to the full pivot grid.
        """
        if not self.valid:
            self.pivot_freq = np.zeros_like(self.pivots)
            return

        self.fe.learn()

        self.pivot_freq = np.zeros_like(self.pivots)
        self.pivot_freq[self.good_pivots] = self.fe.pivot_freq

        # CHANGED: no forced 0/1 outside window — let higher-level logic decide
        # Previously there was special-case clamping by lineage name.


class nested_frequencies(object):
    """
    estimates frequencies of mutually exclusive events such as mutations
    at a particular genomic position or subclades in a tree
    """

    def __init__(self, tps, obs, pivots, **kwargs):
        """
        Parameters
        ----------
        tps : np.array
            array of numerical dates (one per sample in this parent node)
        obs : dict[str, np.array(bool)]
            mapping: key -> boolean array (same length as tps)
        pivots : np.array
            pivot values
        """
        super(nested_frequencies, self).__init__()
        self.tps = np.asarray(tps)
        self.obs = obs
        self.pivots = np.asarray(pivots)
        self.kwargs = kwargs

    def calc_freqs(self):
        """
        Estimate mutually exclusive child frequencies on a shared pivot grid.

        Returns
        -------
        dict
            Mapping of observation key to estimated pivot frequencies.
        """
        # sort by total number of positives per clade (largest first)
        sorted_obs = sorted(self.obs.items(), key=lambda x: x[1].sum(), reverse=True)

        self.remaining_freq = np.ones_like(self.pivots)
        self.frequencies = {}

        # CHANGED: no valid_tps shrinking; use the full tps for all clades
        # This keeps observations aligned and avoids subtle mask bugs.
        for mut, obs in sorted_obs:
            fe = freq_est_clipped(self.tps, obs, self.pivots, name=mut, **self.kwargs)

            if not getattr(fe, "valid", False):
                # no usable data for this clade
                self.frequencies[mut] = np.zeros_like(self.remaining_freq)
                continue

            fe.learn()
            self.frequencies[mut] = self.remaining_freq * fe.pivot_freq
            self.remaining_freq *= (1.0 - fe.pivot_freq)

        return self.frequencies


class tree_frequencies(object):
    """
    Frequency estimator on an abstract lineage tree.

    Assumptions:
    -----------
    - There is a SINGLE global list of observation times stored at the root:
        tree.root.traits['tps'] -> np.array of shape (N,)
      where N is the number of sequences.
    - Each node (internal or leaf) has:
        node.traits['leafs'] -> np.array of integer indices into root.traits['tps']
      corresponding to all samples descended from that lineage.
    - The tree itself is small (lineage-level), but N can be large.
    """

    def __init__(self, tree, timepoints, verbose=0, pc=1e-4, **kwargs):
        """
        Initialize a tree-wide clade frequency estimator.

        Parameters
        ----------
        tree : :class:`baltic.tree.Tree`
            Tree whose nodes carry observation metadata.

        timepoints : int or array-like
            Number of pivots or explicit pivot locations used for frequency
            estimation.

        verbose : int, optional
            Verbosity flag for diagnostic output.

        pc : float, optional
            Pseudocount used to clip frequencies away from 0 and 1.

        **kwargs
            Additional keyword arguments forwarded to nested estimators.
        """

        self.tree = tree
        self.timepoints = timepoints
        self.verbose = verbose
        self.kwargs = kwargs
        self.pc = pc

        self.prepare()

    def prepare(self):
        """
        Prepare:
            - grab global tps from root
            - initialize root frequency trajectory
            - precompute counts per pivot
        """
        # CHANGED: operate on abstract lineage tree with global tps
        indexObs(self.tree.root) ## assign indices to observations

        propagateParentObs(self.tree.root, 'tps') ## propagate observation dates to parents
        propagateParentObs(self.tree.root, 'leafs') ## propagate observation indices to parents

        minT, maxT = min(self.tree.root.traits['tps']), max(self.tree.root.traits['tps'])
        self.pivots = np.linspace(minT, maxT, self.timepoints)

        try:
            self.tps = np.asarray(self.tree.root.traits['tps'])
        except KeyError:
            raise KeyError("Expected tree.root.traits['tps'] with global observation times.")

        if np.isscalar(self.pivots):
            self.frequencies = {self.tree.root.index: np.ones(int(self.pivots))}
        else:
            self.pivots = np.asarray(self.pivots)
            self.frequencies = {self.tree.root.index: np.ones_like(self.pivots)}

        self.counts = count_observations(self.pivots, self.tps)





    def estimate_clade_frequencies(self):
        """
        Compute frequencies for all clades in a top-down manner.

        For each internal node:
          - restrict to its descendant sample indices
          - build a per-child boolean obs array w.r.t. those indices
          - run nested_frequencies
        """

        # internal nodes in pre-order (root first)
        internal_nodes = self.tree.traverse_tree(includeCondition=lambda k: k.is_node())

        for node in internal_nodes:
            if not hasattr(node, 'traits') or 'leafs' not in node.traits:
                continue  # nothing to do

            node_indices = np.asarray(node.traits['leafs'], dtype=int)
            if node_indices.size == 0:
                continue

            node_tps = self.tps[node_indices]

            # build observation dict over children
            obs_to_estimate = {}
            for child in node.children:
                if 'leafs' not in child.traits:
                    continue
                child_indices = np.asarray(child.traits['leafs'], dtype=int)
                if child_indices.size == 0:
                    continue

                # Boolean array per sample in this node: True if sample in this child
                obs = np.isin(node_indices, child_indices)
                if obs.sum() == 0:
                    continue
                obs_to_estimate[child.index] = obs

            if not obs_to_estimate:
                continue

            if self.verbose:
                print(f"Estimating nested frequencies for node {node.index} "
                      f"with {len(obs_to_estimate)} children.")

            ne = nested_frequencies(node_tps, obs_to_estimate, self.pivots,
                                    pc=self.pc, **self.kwargs)
            child_freqs = ne.calc_freqs()

            parent_key = node.index
            if parent_key not in self.frequencies:
                # e.g. root: ensure it has a trajectory (should have from prepare)
                self.frequencies[parent_key] = np.ones_like(self.pivots)

            for clade_key, cfreq in child_freqs.items():
                self.frequencies[clade_key] = self.frequencies[parent_key] * cfreq

        # Assign zero frequencies to tips that did not pass the node_filter.
        # for tip in self.tree.get_external():
        #     if not self.node_filter(tip):
        #         self.frequencies[tip.index] = np.zeros(len(self.pivots))

        return self.frequencies

    def calc_confidence(self):
        """
        Approximate confidence using Bernoulli sampling variance.
        """
        self.confidence = {}
        for key, freq in self.frequencies.items():
            freq = np.asarray(freq)

            self.confidence[key] = np.sqrt(
                (1.0 / (1 + self.counts) + freq * (1 - freq)) /
                (1.0 + self.counts)
            )
        return self.confidence

def clip_freq(ys,threshold=0.00):
    """
    Clip a list of frequencies to the first non-zero and the last non-zero values.
    """
    if sum(ys) == 0.0:
        return []

    xs = range(len(ys))

    # Find first pivot above threshold
    start = next((i for i, y in enumerate(ys) if y >= threshold), None)

    # Find last pivot above threshold
    end_rev = next((i for i, y in enumerate(ys[::-1]) if y >= threshold), None)

    if start is None:
        return []  # nothing above threshold at all

    if end_rev is None:
        # No trailing value above threshold — only leading ones
        end = start
    else:
        end = len(ys) - end_rev - 1  # convert reverse index to forward index

    # ensure at least one pivot
    if start > end:
        end = start

    # avoid slicing at 0-index bug
    if start == 0:
        start = 1

    return list(xs)[start-1:end+1]


def plot_Muller(ax, node, timeline, frequenciesDict, bottom = None, colourFxn = None, outlineColourFxn = None, 
    labelFxn = None, Muller = True, normaliseFreqFxn = None, orientation = 'horizontal', 
    filterFxn = None, clipThreshold = 0.001, level = None, **kwargs):
    """
    Plot a Muller plot on ax starting from a node and accessing its frequency using the size function along a timeline grid.
    """
    from collections import Counter

    valid_orientations = ['vertical', 'horizontal']
    assert orientation in valid_orientations, f"Orientation '{orientation}' not recognised. Must be one of {', '.join(valid_orientations)}."
    
    if level is None: level=1
    if filterFxn is None: filterFxn=lambda n: True

    if colourFxn is None: colourFxn = lambda k: 'lightgray'
    if outlineColourFxn is None: outlineColourFxn = lambda k: 'white'
    
    localKwargs = dict(kwargs)
    
    if 'facecolor' in localKwargs or 'fc' in localKwargs:
        warnings.warn("facecolor/fc provided in kwargs will be ignored. Muller plot colours are controlled via colourFxn.")
    if 'edgecolor' in localKwargs or 'ec' in localKwargs:
        warnings.warn("edgecolor/ec provided in kwargs will be ignored. Muller plot outlines are controlled via outlineColourFxn.")
    if 'zorder' in localKwargs:
        warnings.warn("zorder provided in kwargs will be ignored. Muller plot needs to use zorder to ensure descendant lineages are displayed on top.")

    ys=frequenciesDict[node.index] ## returns values of a node over time
    if normaliseFreqFxn: ## if normalise function available - normalise frequencies
        ys=normaliseFreqFxn(ys)
    
    bottom=[0-y/2 for y in ys] if bottom==None else bottom ## if no bottom provided - compute it 
    
    xs=range(len(ys)) ## indices of frequencies
    xx=clip_freq(ys, clipThreshold) ## reindex frequencies to clip zeroes at the beginning and end
    
    localKwargs['fc'] = colourFxn(node)
    localKwargs['ec'] = outlineColourFxn(node)
    localKwargs['zorder'] = level
    if 'alpha' not in localKwargs: localKwargs['alpha'] = 1.0
    
    clipped_timeline = [timeline[t] for t in xx]
    clipped_bottom = [bottom[t] for t in xx]
    clipped_values = [bottom[t] + ys[t] for t in xx]

    label = labelFxn(node) if labelFxn else ''

    
    if orientation == 'horizontal':
        fillBetweenFxn = ax.fill_between
    elif orientation == 'vertical':
        fillBetweenFxn = ax.fill_betweenx

    if filterFxn(node):
        fillBetweenFxn(clipped_timeline, clipped_bottom, clipped_values, label = label, **localKwargs) ## plot frequency
    
    if node.is_node(): ## if node has children
        children = node.children ## children of node
        
        clipped_children = [list(clip_freq(frequenciesDict[ch.index], clipThreshold)) for ch in children] ## clip frequencies of child
        clipped_freqs = Counter(sum(clipped_children, [])) ## count how many timeline indices are left
        
        children_sum=[sum([frequenciesDict[ch.index][t] for ch in children]) for t in xs] ## total frequency of children
        
        if normaliseFreqFxn:
            children_sum = normaliseFreqFxn(children_sum) ## if normalise available - normalise child frequencies
        
        available_space = [(ys[t] - children_sum[t]) for t in xs] ## node frequency - all children frequencies = padding space left
        N_children = [clipped_freqs[t] + 1 if clipped_freqs[t] > 0 else 1 for t in xs] ## count+1 of indices = number of children at any point
        unique_space = [available_space[t] / N_children[t] for t in xs] ## how much padding to add
        
        temp_bottom = bottom ## start with bottom
    
        for c,child in enumerate(children): ## iterate over children
            
            child_freq = clip_freq(frequenciesDict[child.index], clipThreshold) ## clip child frequencies

            if Muller:
                padded_bottom = [temp_bottom[t] + unique_space[t] if t in child_freq else temp_bottom[t] for t in xs] ## pad bottom with space available, if a child is present
            else:
                padded_bottom = [temp_bottom[t] if t in child_freq else temp_bottom[t] for t in xs]
            
            
            temp_bottom = plot_Muller(ax = ax, node = child, timeline = timeline, 
                               frequenciesDict = frequenciesDict, bottom = padded_bottom, 
                               colourFxn = colourFxn, outlineColourFxn = outlineColourFxn, 
                               labelFxn = labelFxn, Muller = Muller, normaliseFreqFxn = normaliseFreqFxn, 
                               filterFxn = filterFxn, orientation = orientation, clipThreshold = clipThreshold, **kwargs) ## draw frequency for each child, padding as you go
        
    return [bottom[t] + ys[t] for t in xs] ## new bottom is bottom+this node's values


#### convex hulls around clades
#### matrix + tree
#### LBI
#### layering function for colouring (e.g. marking discovery paths)


def position_Bezier_control(pointA, pointB, height, frac):
    """ 
    Given a line defined by 2 points A & B, 
    find a third point at a given distance (height) that defines a line perpendicular to line AB which intercepts AB at fraction (frac) along AB.
    Equation derived by Luiz Max Fagundes de Carvalho (University of Edinburgh).
    """
    x1, y1 = pointA
    x2, y2 = pointB

    sign = 1
    if x1 > x2:
        sign = -1

    slope = (y2 - y1) / (x2 - x1)
    d=np.sqrt((y2 - y1)**2 + (x2 - x1)**2) ## distance between points
    
    h=np.sqrt(height**2 + (d * frac)**2) ## distance between desired height and point along line

    n1 = x1 + h * np.cos(np.arctan(height / float(d) / frac) + np.arctan(slope)) * sign ## magic
    n2 = y1 + h * np.sin(np.arctan(height / float(d) / frac) + np.arctan(slope)) * sign

    return (n1, n2) ## return third point's coordinate
##########################
def plot_root_to_tip(ax, tree,
                     colour = None, colourFxn = None,
                     pointSize = None, pointSizeFxn = None,
                     outline = True, outlineColour = None, outlineColourFxn = None,
                     outlineSize = None, outlineSizeFxn = None,
                     orientation = 'horizontal',
                     targetFxn = None,
                     plotRegression = True, plotTree = False,
                     highlightOutliers = False, outlierThres = None, outlierColour = None,
                     pointKwargs = {}, lineKwargs = {}, treeKwargs = {}):

    import math
    from baltic.bt_utils import desaturate_cmap, calendar_to_decimal_date, decimal_to_calendar_date
    from scipy.stats import linregress

    assert tree.treeType == 'divergence', f"Cannot run root-to-tip regression on a time tree. Branch lengths must be in substitution space."
    assert orientation in ['horizontal', 'vertical'], f"Orientation {orientation} not recognised. Must be 'horizontal' or 'vertical'."

    func_logger = logging.getLogger("baltic.bt_utils.plot_root_to_tip")
    func_logger.setLevel(logging.INFO)
    func_logger.propagate = False
    if not func_logger.handlers:
        import sys
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(logging.Formatter("[root_to_tip] %(message)s"))
        handler.setLevel(logging.INFO)
        func_logger.addHandler(handler)

    if targetFxn is None: targetFxn = lambda k: k.is_leaf()

    localPointKwargs = dict(pointKwargs)
    localLineKwargs = dict(lineKwargs)
    localTreeKwargs = dict(treeKwargs)

    tipOrder = tree.get_external(targetFxn)

    ys = [k.height for k in tipOrder]
    xs = [k.absoluteTime for k in tipOrder if k.absoluteTime]

    assert len(xs) > 0, f"Tips do not have an assigned `absoluteTime` attribute that represents tip collection dates."
    assert len(xs) == len(ys), f"Some branches missing collection dates. Branhes with heights: {len(ys)}; with collection dates: {len(xs)}."
    #################

    # Colour logic
    if colour is not None and colourFxn is not None:
        raise ValueError(
            "Cannot specify both colour and colourFxn. Please use only one."
        ) ## should be a warning, since this eventuality is handled in the next line
    if colour is None and colourFxn is None:
        cmap = desaturate_cmap(mpl.cm.Spectral, 0.6)
        norm = mpl.colors.Normalize(min(xs), max(xs))
        colourFxn = lambda k: cmap(norm(k.absoluteTime))

    elif colourFxn is None:
        colourFxn = lambda k: colour

    if outlineColour is None and outlineColourFxn is None:
        outlineColourFxn = lambda k: "k"
    elif outlineColourFxn is None:
        outlineColourFxn = lambda k: outlineColour

    # Point size logic
    if pointSize is not None and pointSizeFxn is not None:
        raise ValueError(
            "Cannot specify both pointSize and pointSizeFxn. Please use only one."
        )
    if pointSize is None and pointSizeFxn is None:
        pointSizeFxn = lambda k: 60
    elif pointSizeFxn is None:
        pointSizeFxn = lambda k: pointSize

    # Outline size logic
    if outlineSize is not None and outlineSizeFxn is not None:
        raise ValueError(
            "Cannot specify both outlineSize and outlineSizeFxn. Please use only one."
        )
    if outlineSize is None and outlineSizeFxn is None:
        outlineSizeFxn = lambda k: ((2+math.sqrt(pointSizeFxn(k)/math.pi))**2)*math.pi
    elif outlineSizeFxn is None:
        outlineSizeFxn = lambda k: outlineSize

    slope, intercept, rval, _, std_err = linregress(xs, ys) ## run linear regression

    tmrca = -intercept / slope
    func_logger.info(f"Root-to-tip regression parameters: slope = {slope:.4e} intercept = {intercept:.4g} TMRCA = {tmrca:.3f} "
                f"({decimal_to_calendar_date(tmrca, fmt='%Y-%b-%d')})")

    sizes = [pointSizeFxn(k) for k in tipOrder]
    cs = [colourFxn(k) for k in tipOrder]

    if outlierThres is None: ## default outlier threshold is residual standard deviation
        residuals = {k.name: float(k.height - (k.absoluteTime * slope + intercept)) for k in tipOrder}
        justResiduals = [residuals[k.name] for k in tipOrder]
        sigma = np.sqrt(np.sum(np.pow(justResiduals,2)) / (len(tipOrder) - 2))
        outlierThres = 2 * sigma

    outliers = [k for k in tipOrder if abs(k.height - (k.absoluteTime * slope + intercept)) >= outlierThres]

    if plotRegression:

        x_grid = np.linspace(min(xs), max(xs))
        y_grid = [x * slope + intercept for x in x_grid]

        if orientation == 'vertical':
            x_grid, y_grid = y_grid, x_grid

        if 'color' not in localLineKwargs: localLineKwargs['color'] = 'k'
        if 'lw' not in localLineKwargs and 'linewidth' not in localLineKwargs: localLineKwargs['lw'] = 3
        if 'ls' not in localLineKwargs and 'linestyle' not in localLineKwargs: localLineKwargs['ls'] = '--'
        ax.plot(x_grid, y_grid, **localLineKwargs)

    if orientation == 'vertical':
        xs, ys = ys, xs

    if 'zorder' not in localPointKwargs: localPointKwargs['zorder'] = 4
    if 'edgecolor' not in localPointKwargs: localPointKwargs['edgecolor'] = 'none'
    ax.scatter(
            xs,
            ys,
            s=sizes,
            facecolor=cs,
            **localPointKwargs,
        )  ## put a circle at each tip
    if outline:
        os = [outlineSizeFxn(k) for k in tipOrder]
        outlierColour = 'darkred' if outlierColour is None else outlierColour
        ocs = [outlierColour if (highlightOutliers and k in outliers) else outlineColourFxn(k) for k in tipOrder]
        localPointKwargs['zorder'] -= 1 ## decrease zorder for outlines
        ax.scatter(
            xs,
            ys,
            s=os,
            facecolor=ocs,
            **localPointKwargs,
        )  ## put a circle at each tip

    for k in tipOrder:
        if np.diff(k.absoluteTimeRange)[0] > 1e-10:
            uxs, uys = k.absoluteTimeRange, [k.height, k.height]
            if orientation == 'vertical':
                uxs, uys = uys, uxs
            ax.plot(uxs, uys, color = 'k')

    if plotTree:

        xFxn = lambda k: k.absoluteTime if k.is_leaf() else (k.height - intercept) / slope
        yFxn = lambda k: k.height

        if orientation == 'vertical':
            xFxn, yFxn = yFxn, xFxn

        if 'connectionType' in localTreeKwargs and localTreeKwargs['connectionType'] != 'direct':
            warnings.warn(f"plot_root_to_tip does not support connectionType other than 'direct'.")
        localTreeKwargs['connectionType'] = 'direct'
        if 'alpha' not in localTreeKwargs: localTreeKwargs['alpha'] = 0.05

        tree.plot_tree(ax, xCoordinateFxn = xFxn, yCoordinateFxn = yFxn, **localTreeKwargs)

    return (outliers, slope, intercept, residuals)



def plot_skygrid(ax, logFile, burnin=None, mostRecent=None, hpdLvl=0.95, orientation='horizontal', logAxis=True, plotRootHPD=True, **kwargs):
    localKwargs = dict(kwargs)

    import csv
    from baltic.bt_utils import hpd
    
    assert mostRecent is None or isinstance(mostRecent, float) or isinstance(mostRecent, int), f"mostRecent value '{mostRecent}' not recognised. Must be None, float or int"
    assert orientation in ['horizontal', 'vertical'], f"Orientation {orientation} not recognised. Must be 'horizontal' or 'vertical'"

    func_logger = logging.getLogger("baltic.bt_utils.plot_skygrid")
    func_logger.setLevel(logging.INFO)
    func_logger.propagate = False

    if not func_logger.handlers:
        import sys
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(logging.Formatter("[plot_skygrid] %(message)s"))
        handler.setLevel(logging.INFO)
        func_logger.addHandler(handler)

    if burnin == None:
        burnin = 10e6
        # logger.warning('No burnin set, defaulting to 10M states.')

    handle = open(logFile,'r')
    reader = csv.DictReader((line for line in handle if line.startswith('#') == False), delimiter = '\t')
    header = reader.fieldnames

    assert any(k.startswith('skygrid.logPopSize') for k in header), "Missing skygrid.logPopSize columns"
    assert 'skygrid.cutOff' in header, "Missing skygrid.cutOff column"

    def compute_skygrid_interval(lineDict):
        """
        Computes root date and most recent date from log file or provided most recent date (otherwise time is years before most recent date)
        """
        rootHeightCandidates = ['treeModel.rootHeight', 'rootHeight']

        rootDate, mostRecentDate = None, None
        
        if mostRecent is None and 'age(root)' in lineDict: ## determining most recent date from rootDate
            rootDate = float(lineDict['age(root)'])
            mostRecentDate = rootDate + next((float(lineDict[candidate]) for candidate in rootHeightCandidates if candidate in lineDict), None)
            assert mostRecentDate is not None, f"Unable to identify most recent tip date from log file provided. None of the expected parameters required for calculation ({', '.join(rootHeightCandidates)}) found in log file."
            
        elif mostRecent is None or isinstance(mostRecent, float) or isinstance(mostRecent, int): ## 
            rootDate = next((float(lineDict[candidate]) for candidate in rootHeightCandidates if candidate in lineDict), None)
            assert rootDate is not None and rootDate != '', f"Unable to identify root / tree height from log file provided. None of the expected parameters {', '.join(rootHeightCandidates)} found in log file."
            
            if mostRecent is None:
                rootDate *= -1 ## not most recent date provided, root date will be in years before most recent tip
                mostRecentDate = 0.0
            else:
                mostRecentDate = mostRecent
        
        return rootDate, mostRecentDate
    
    skygridPopSizes = sorted([key for key in header if key.startswith('skygrid.logPopSize')], key = lambda f: int(f.replace('skygrid.logPopSize', '')))

    skygrid = {key: [] for key in skygridPopSizes}
    rootAge = []
    cutoff = None

    storeMostRecent = None
    
    for l in reader:
        state = int(l['state'])
        if state >= burnin:
            for key in skygridPopSizes:
                popSize = float(l[key])
                skygrid[key].append(popSize)

            rootAgeParam, mostRecentDate = compute_skygrid_interval(l)

            if storeMostRecent:
                assert mostRecentDate == storeMostRecent, f"Estimated recent date ({mostRecentDate}) differs from previous MCMC state ({storeMostRecent}). This can happen when tip dates sampled during MCMC become the most recent tip dates. Check the XML that produced the log file."
            rootAge.append(rootAgeParam)

            if cutoff is None:
                cutoff = float(l['skygrid.cutOff'])

            storeMostRecent = mostRecentDate
    handle.close()
    
    start = mostRecentDate - cutoff
    end = mostRecentDate

    xs = np.linspace(start, end, len(skygridPopSizes))

    ys = [np.pow(np.e, np.mean(skygrid[key])) for key in skygridPopSizes[::-1]]
    y_lower, y_upper = zip(*[np.pow(np.e, hpd(skygrid[key], hpdLvl)) for key in skygridPopSizes[::-1]])

    fillPlotFxn = ax.fill_between
    if orientation == 'vertical':
        fillPlotFxn = ax.fill_betweenx

    if 'facecolor' not in localKwargs and 'fc' not in localKwargs: localKwargs['facecolor'] = 'lightskyblue'
    fillPlotFxn(xs, y_lower, y_upper, zorder = 9, **localKwargs)
    if orientation == 'vertical':
        xs, ys = ys, xs

    ax.plot(xs, ys, color = 'k', lw = 2, zorder = 10) ## plot mean

    if plotRootHPD: ## plotting root HPD lines
        rootPlotFxn = ax.axvline
        if orientation == 'vertical':
            rootPlotFxn = ax.axhline

        rootMean = np.mean(rootAge)
        rootLo, rootHi = hpd(rootAge, hpdLvl)

        if isinstance(mostRecent, float) or isinstance(mostRecent, int):
            rootMean = end - rootMean
            rootLo = end - rootLo
            rootHi = end - rootHi

        rootPlotFxn(rootMean, color = 'dimgray', zorder = 10)
        rootPlotFxn(rootLo, color = 'dimgray', ls = '--', zorder = 10)
        rootPlotFxn(rootHi, color = 'dimgray', ls = '--', zorder = 10)

    if logAxis:
        if orientation == 'horizontal':
            ax.set_yscale('log')
        elif orientation == 'vertical':
            ax.set_xscale('log')

    return ax


def connect_tree_to_map(treeAx, mapAx, tree, tipCoordinates, destinationProjection, targetFxn=None, shoulderPositionFxn=None, colourFxn=None, originProjection=None, **lineKwargs):
    """
    treeAx is the Axes object where the tree is plotted already.
    mapAx is the Axes object where the map is plotted already.
    tree is a baltic Tree object.
    tipCoordinates is a dict of {tip_name: (x, y)} coordinates.
    destinationProjection is a cartopy ccrs projection object used in mapAx Axes.
    targetFxn is a filtering function for tree tips (default: always True).
    shoulderPositionFxn is a function that connects the x-axis coordinate provided by this function to the map (default: None, tips are directly connected to the map).
    colourFxn returns the colour of a tip.
    originProjection is a cartopy ccrs object used to transform the coordinate system of tipCoordinates (default: ccrs.PlateCarree()).
    """
    from matplotlib.patches import ConnectionPatch
    import cartopy.crs as ccrs
    from baltic.bt_utils import desaturate_cmap
    
    if originProjection is None: originProjection = ccrs.PlateCarree()
    if targetFxn is None: targetFxn = lambda k: True
    
    cmap = desaturate_cmap(mpl.cm.Spectral, 0.6)
    if colourFxn is None: colourFxn = lambda k: cmap(k.y / tree.ySpan)

    localLineKwargs = dict(lineKwargs)

    
    for k in tree.get_external(targetFxn):
        assert k.name in tipCoordinates, f"Tip {k.name} coordinates not found in tipCoordinates dict"

        fc = colourFxn(k)
        
        lat, lon = tipCoordinates[k.name]
        latT, lonT = destinationProjection.transform_point(lon, lat, src_crs=originProjection) ## convert to plotted map coordinates

        if shoulderPositionFxn is None:
            x, y = k.x, k.y
        else:
            x, y = shoulderPositionFxn(k), k.y
        
        connection = ConnectionPatch(xyA=(x, y), 
                                     coordsA=treeAx.transData, 
                                     axesA=treeAx, 
                                     xyB=(latT, lonT), 
                                     coordsB=mapAx.transData, 
                                     axesB=mapAx, 
                                     color=fc, **localLineKwargs) ## colour, line style, linewidth, order, transparency
        #### check if connection will be present
        axA = connection.axesA
        axB = connection.axesB
        pA = axA.transData.transform(connection.xy1)
        pB = axB.transData.transform(connection.xy2)
        bbA = axA.bbox
        bbB = axB.bbox
    
        if (bbA.contains(*pA) or bbB.contains(*pB)) and shoulderPositionFxn is not None: ## checks if line will be visible
            treeAx.plot([k.x, x], [y, y], color=fc, **localLineKwargs) ## straight line from the tip to the x-axis coordinate provided by the shoulderPositionFxn
        
        mapAx.add_patch(connection)

    return treeAx, mapAx

def plot_tangled_chain(ax, treeList, colourDict=None, padding=None, treeSpaceFxn=None, treeSpace=None, normaliseY=True, treeKwargs={}, pointKwargs={}, **kwargs):
    from matplotlib.collections import LineCollection

    localKwargs = dict(kwargs)
    localTreeKwargs = dict(treeKwargs)
    localPointKwargs = dict(pointKwargs)

    if treeSpace is not None and treeSpaceFxn is not None: ## treeSpace is how much space is left between consecutive trees, takes current tree if specified as a function
        raise ValueError(
            "Cannot specify both treeSpace and treeSpaceFxn. Please use only one."
        ) ## should be a warning, since this eventuality is handled in the next line
    if treeSpace is None and treeSpaceFxn is None:
        treeSpaceFxn = lambda k: treeList[0].treeHeight * 0.20 ## 20% of first tree height is space between all trees
    elif treeSpaceFxn is None:
        treeSpaceFxn = lambda k: treeSpace

    if padding is None: ## padding is proportion of treeSpace protrudes beyond previous tree and before next tree (0 == line finishes at last tip of current tree and goes to root of next, 0.5 == line goes to middle between consecutive trees and switches abruptly)
        padding = 0.1
    else:
        assert 0.0 <= padding <= 0.5, f"Padding (given as {padding}) should be a float between 0 and 0.5."

    if colourDict is None: ## colourMap is dict that assigns colours to tips according to their y-axis order in first tree

        colourDict = {}

        cmap = desaturate_cmap(mpl.cm.Spectral, 0.6)
        firstTreeTips = treeList[0].get_external()

        for i,k in enumerate(sorted(firstTreeTips, key = lambda q: q.y)):
            colourDict[k.name] = cmap(i / (len(firstTreeTips) - 1))
    else:
        assert isinstance(colourDict, dict), f"colourDict must be class dict, not {type(colourDict)}"
    
    if 'coordinateFxn' in localTreeKwargs: logger.warning(f"Custom x coordinate function for tree was specified but will be overriden for tangled chain visualisation.")
    # if 'xCoordinateFxn' not in localTreeKwargs: localTreeKwargs['xCoordinateFxn'] = lambda k: k.x + cumulativeX

    if len(localPointKwargs) > 0: ## tip points are required - override xCoordinateFxn, assign default colours if nothing specified
        # if 'xCoordinateFxn' in localPointKwargs: logger.warning(f"Custom x coordinate function for points was specified but will be overriden for tangled chain visualisation.")
        # localPointKwargs['xCoordinateFxn'] = lambda k: k.x + cumulativeX
        
        if 'colour' not in localPointKwargs and 'colourFxn' not in localPointKwargs:
            logger.warning(f"Point colours were not specified, defaulting to tangled chain colour defaults. This may cause issues if targetFxn is not set to identify tips.")
            localPointKwargs['colourFxn'] = lambda k: colourMap[k.name]

    ## compute tree coordinates in tangled chain organisation
    cumulativeX = 0 ## tracks x coordinate as we plot consecutive trees
    for t in range(len(treeList)):
        spaceUnit = treeSpaceFxn(treeList[t])

        treeList[t]._assign_tree_coordinates()

        for k in treeList[t].Objects:
            k.x = cumulativeX + k.x

        cumulativeX += treeList[t].treeHeight + spaceUnit ## increment x-axis


    connectionCoordinates = []
    connectionColours = []

    xCoordinateFxn = lambda k: k.x
    cumulativeX = 0

    for curTree, nexTree in zip(treeList, treeList[1:]): ## iterate over pairs of consecutive trees    

        if normaliseY:
            yCoordinateFxn = lambda k: k.y / curTree.ySpan
        else:
            yCoordinateFxn = lambda k: k.y

        curTree.plot_tree(ax, recomputeCoordinates=False, autoSort=False, xCoordinateFxn=xCoordinateFxn, yCoordinateFxn=yCoordinateFxn, **localTreeKwargs) ## plot current tree

        if len(localPointKwargs) > 0:
            curTree.plot_points(ax, recomputeCoordinates=False, xCoordinateFxn=xCoordinateFxn, yCoordinateFxn=yCoordinateFxn, **localPointKwargs) ## add points if specified

        for curTip in curTree.get_external(): ## iterate over tips in current tree
            c = colourDict[curTip.name] if curTip.name in colourDict else 'lightgray'

            nexTip = nexTree.get_external(filterFxn = lambda k: k.name == curTip.name) ## identify matching tip
            if len(nexTip)>0:
                nexTip = nexTip[0]
                curX = curTip.x
                curY = curTip.y

                lineAfterX = cumulativeX + curTree.treeHeight + spaceUnit * padding
                lineBeforeX = cumulativeX + curTree.treeHeight + spaceUnit * (1 - padding)

                nexX = nexTip.x
                nexY = nexTip.y

                if normaliseY: ## tip position is fraction between 0 and 1
                    curY /= curTree.ySpan
                    nexY /= nexTree.ySpan

                connectionCoordinates.append([(curX, curY),
                                              (lineAfterX, curY),
                                              (lineBeforeX, nexY),
                                              (nexX, nexY)]) ## coordinates of tangled line
                connectionColours.append(c) ## colour of tangled line

        cumulativeX += curTree.treeHeight + spaceUnit

    if normaliseY:
        yCoordinateFxn = lambda k: k.y / nexTree.ySpan

    nexTree.plot_tree(ax, recomputeCoordinates=False, autoSort=False, xCoordinateFxn=xCoordinateFxn, yCoordinateFxn=yCoordinateFxn, **localTreeKwargs) ## plot last tree
    if len(localPointKwargs) > 0:
        nexTree.plot_points(ax, recomputeCoordinates=False, xCoordinateFxn=xCoordinateFxn, yCoordinateFxn=yCoordinateFxn, **localPointKwargs) ## plot its points

    if 'zorder' not in localKwargs: localKwargs['zorder'] = 0
    ax.add_collection(LineCollection(connectionCoordinates, color=connectionColours, **localKwargs)) ## add tangled lines

    return treeList

def plot_gradient_clade_tree(ax, tree, designatedNodes=None, nodeDesignationFxn=None,
    colour=None, colourFxn=None,
    controlPointFraction=0.1,
    tipLen=0.5,
    padY=0.5,
    outlineColour=None, outlineColourFxn=None,
    minAlpha=0.0, maxAlpha=0.4,
    outline=True,
    precision=50,
    outlineKwargs=None, treeKwargs=None, cladeBranchKwargs=None,
    ):
    """
    Plot tree into axes with designated nodes (via list or function) being shown as a gradient.

    Inspired by JT McCrone's visualisation: https://www.nature.com/articles/s41586-022-05200-3
    Initially implemented by Karthik Gangavarapu.
    """
    from matplotlib.patches import Polygon
    from matplotlib.collections import LineCollection
    from baltic.bt_utils import draw_gradient_polygon, five_point_bezier

    if designatedNodes is not None and nodeDesignationFxn is not None:
        raise ValueError(
            "Cannot specify both designatedNodes and nodeDesignationFxn. Please use only one."
        )
    if designatedNodes is None and nodeDesignationFxn is None:
        raise ValueError(
            "Must specify either designatedNodes or nodeDesignationFxn."
        )
    if nodeDesignationFxn is None:
        nodeDesignationFxn = lambda k: k in designatedNodes
    
    if colour is not None and colourFxn is not None:
        raise ValueError(
            "Cannot specify both colour and colourFxn. Please use only one."
        ) ## should be a warning, since this eventuality is handled in the next line
    if colour is None and colourFxn is None:
        colourFxn = lambda k: "k"
    elif colourFxn is None:
        colourFxn = lambda k: colour

    if outlineColour is not None and outlineColourFxn is not None:
        raise ValueError(
            "Cannot specify both outlineColour and outlineColourFxn. Please use only one."
        ) ## should be a warning, since this eventuality is handled in the next line
    if outlineColour is None and outlineColourFxn is None:
        outlineColourFxn = lambda k: colourFxn(k)
    elif outlineColourFxn is None:
        outlineColourFxn = lambda k: outlineColour
    
    assert 0.0 <= maxAlpha <= 1.0 and 0.0 <= minAlpha <= 1.0, f"minAlpha ({minAlpha}) and maxAlpha ({maxAlpha}) must be in range [0.0, 1.0]."
    assert maxAlpha >= minAlpha, f"maxAlpha must be greater than minAlpha. Currently minAlpha is {minAlpha} and maxAlpha is {maxAlpha}"
    
    localOutlineKwargs = dict(outlineKwargs) if outlineKwargs else dict()
    localTreeKwargs = dict(treeKwargs) if treeKwargs else dict()
    localCladeBranchKwargs = dict(cladeBranchKwargs) if cladeBranchKwargs else dict()
    if 'linewidth' not in localCladeBranchKwargs and 'lw' not in localCladeBranchKwargs: localCladeBranchKwargs['lw'] = 2
    assert 'color' not in localCladeBranchKwargs, f"Not allowed to specify 'color' in cladeBranchKwargs, branch colour in gradient clades is handled by either colour or colourFxn parameters."
    
    gradientClades = tree.get_internal(nodeDesignationFxn)
    assert len(gradientClades) > 0, f"No nodes specified."
    
    gradientSubtrees = set()
    done = set()
    
    for node in sorted(gradientClades, key=lambda k: len(k.leaves)): ## iterate over designated nodes starting from smallest clades

        startX = node.absoluteTime
        startY = node.y
        
        subtree = tree.traverse_tree(node, includeCondition=lambda k: True)
        gradientSubtrees = gradientSubtrees.union(subtree)
        descendants = [leaf for leaf in subtree if leaf.is_leaflike()] ## only keep leaves
        
        if len(descendants) == 0:
            return ax
        
        yLo, yHi = node.yRange
        yLo -= padY
        yHi += padY
    
        endX = max([k.absoluteTime for k in descendants]) ## endX
        
        cladeLen = endX - startX

        startingPoint = (startX, startY) ## where clade begins
        earlyX = startX + cladeLen * controlPointFraction * 0.5
        midX = startX + cladeLen * controlPointFraction
        
        cpLo = [startingPoint, 
                (earlyX, yLo), 
                (earlyX, yLo), 
                (midX, yLo), 
                (endX, yLo)]
        
        cpHi = [startingPoint, 
                (earlyX, yHi), 
                (earlyX, yHi), 
                (midX, yHi), 
                (endX, yHi)]
        
        loXY = list(zip(*five_point_bezier(cpLo, precision=precision)))
        hiXY = list(zip(*five_point_bezier(cpHi, precision=precision)))
        
        xy = np.vstack([[startX, startY], loXY, hiXY[::-1]]) ## gradient clade coordinates: starting point -> lower fan -> upper fan (reversed coords)

        draw_gradient_polygon(ax=ax, polygonXY=xy, extent=[startX, endX, yLo, yHi], axis='x', reverse=True, 
                              colour=colourFxn(node), minAlpha=minAlpha, maxAlpha=maxAlpha)
    
        if outline: ## plotting outline
            assert 'color' not in localOutlineKwargs, f"Not allowed to specify 'color' in outlineKwargs, outline colour is handled by either outlineColour or outlineColourFxn parameters."
            ax.plot(*list(zip(*loXY)), color=outlineColourFxn(node), **localOutlineKwargs)
            ax.plot(*list(zip(*hiXY)), color=outlineColourFxn(node), **localOutlineKwargs)
        ###############
        branches = [] ## contains branch coordinates
        colours = [] ## contains branch colours

        assert 0.0 <= tipLen <= 1.0, f"Provided tipLen ({tipLen}) must be within interval [0.0, 1.0]."
        fs = np.linspace((1 - tipLen), 1, precision) ## line segments to be plotted with decreasing transparency
        
        for b in descendants: ## iterate over node's descendants
            if b in done: ## only want branches not yet plotted
                continue
            
            x, y = b.absoluteTime, b.y
            parX = node.absoluteTime
            coords = [(x - (x - parX) * (f - (1 - tipLen)), y) for f in fs[::-1]]
    
            z = np.empty((precision - 1, 4))
            rgb = mpl.colors.colorConverter.to_rgb(colourFxn(b))
            z[:,  :3] = rgb
            z[:,-1] = np.append((np.logspace(0, 1, precision - 2) - 1) / 10, 1)
            
            branches.extend([(c1, c2) for c1, c2 in zip(coords, coords[1:])])
            colours.extend(z)
            
            done.add(b) ## designate branch as done
        ax.add_collection(LineCollection(branches, color=colours, zorder=1, **localCladeBranchKwargs))
    
    tree.plot_tree(ax=ax, targetFxn=lambda k: k not in gradientSubtrees or nodeDesignationFxn(k), colourFxn=colourFxn, connectionType='elbow', autoSort=False, zorder=1, **localTreeKwargs)

    return ax

def plot_height_95hpds(ax, tree, targetFxn=None, traitName='height_95%_HPD', width=0.5, lastTipDate=None, **kwargs):
    """
    Plot 95% height HPDs on branches with trait.
    """
    from matplotlib.patches import Rectangle

    localKwargs = dict(kwargs)

    if 'fc' not in localKwargs and 'facecolor' not in localKwargs: localKwargs['facecolor'] = 'lightgray'
    if 'ec' not in localKwargs and 'edgecolor' not in localKwargs: localKwargs['edgecolor'] = 'dimgray'
    if 'zorder' not in localKwargs: localKwargs['zorder'] = 0
    if 'alpha' not in localKwargs: localKwargs['alpha'] = 0.4
    
    if targetFxn is None: targetFxn = lambda k: True
    if lastTipDate is None: lastTipDate = tree.mostRecent
    
    for k in tree.Objects:
        if targetFxn(k) == False or traitName not in k.traits:
            continue

        x, y = k.absoluteTime, k.y
        lo, hi = [lastTipDate - stat for stat in k.traits[traitName]] ## expect two entries
        barWidth = hi - lo
        barHeight = width ## more intuitive to think of bars as having width, but it'll correspond to height parameter of mpl Rectangle
        
        ax.add_patch(Rectangle((lo, y - barHeight / 2), width=barWidth, height=barHeight, **localKwargs))
    
    return ax

def plot_reticulations(ax, tree, excludeFxn=None, colour=None, colourFxn=None, plotSegMatrix=False, segMatrixDist=1.1, segNames=None, segOrder=None, segWidth=1, segHeight=1, **kwargs):
    """
    Add a little matrix at the end of the reassortment network to indicate which segments are reassorting.
    """
    from matplotlib.collections import LineCollection

    if excludeFxn is None: excludeFxn = lambda k: False
    
    if colour is not None and colourFxn is not None:
        raise ValueError(
            "Cannot specify both colour and colourFxn. Please use only one."
        ) ## should be a warning, since this eventuality is handled in the next line
    if colour is None and colourFxn is None:
        colourFxn = lambda k: "gray"
    elif colourFxn is None:
        colourFxn = lambda k: colour
    
    if plotSegMatrix == False and segNames is not None: ## check for superfluous parameters
        logger.warning(f"Provided segment names {segNames} will be ignored because plotSegMatrix is set to False.")
    
    if plotSegMatrix: ## identifies segment numbers in network, checks if names are available for renaming
        import re
        from matplotlib.patches import Rectangle
        
        segFind = re.compile('^seg([0-9]+)$') ## typical CoalRe formatting of segment names
        segments = []
        for trait in tree.root.traits: ## iterate over root's traits
            match = segFind.match(trait) ## check if trait key is formatted like a CoalRe segment
            if match:
                segments.append(int(match.group(1))) ## remember segment

        if segNames is not None: ## check that provided dicts are the right format
            assert isinstance(segNames, dict), f"segNames should be dict ({type(segNames)} provided) with integers [0, 1, ..., N] as keys mapping trait entries 'seg0', 'seg1', ..., 'segN' in reassortment network to segment names."
            assert all(seg in segNames for seg in segments), f"Could not find names for following segments in segNames: {', '.join(str(seg) for seg in segments if seg not in segNames)}"
            if segOrder is None:
                segOrder = sorted(segNames.values())
            else:
                assert len(segOrder) == len(segNames), f"segOrder has different number of entries ({len(segOrder)}) than segNames ({len(segNames)})"
                assert all(seg in segOrder for seg in segNames.values()), f"Could not find order for following segments in segOrder: {', '.join(str(seg) for seg in segNames.values() if seg not in segOrder)}"

    branches = []
    colours = []

    pointerLine = []
    pointerLineColours = []
    
    treeEnd = tree.root.absoluteTime + tree.treeHeight * segMatrixDist ## get x coordinate at some position after the tree

    reassortmentCounter = 0 ## count reassortments
    
    for k in sorted(tree.Objects, key=lambda q: -q.y): ## iterate over all branches from top to bottom along y-axis
        if isinstance(k, Reticulation) and excludeFxn(k) == False: ## is a reassortment edge or excluded
            x, y = k.absoluteTime, k.y ## get reassortment origin coordinate
            xd, yd = k.target.absoluteTime, k.target.y ## get reassortment destination coordinate

            branches.append(((x, y), (xd, yd))) ## add to list of branches

            marker = '^' if yd > y else 'v' ## marker up if reticulation target (destination) is above origin, downward otherwise

            fc = colourFxn(k)
       
            ax.scatter(xd, yd, s=50, fc='w', ec='none', marker=marker, zorder=10) ## add arrow-like ending of reassortment edge
            ax.scatter(xd, yd, s=120, fc='k', ec='none', marker=marker, zorder=9)

            colours.append('k')

            reassortmentCounter += 1
            
            if plotSegMatrix:
                
                pointerLine.append(((x, y), (treeEnd, y)))
                pointerLineColours.append('gray')

                for i, seg in enumerate(segOrder): ## iterate over segments
                    sx = treeEnd + i * segWidth ## find x position, offset by segment index

                    segMatch = [segIdx for segIdx in segments if segNames[segIdx] == seg][-1] ## identify segment in trait dict
                    segInvolved = k.traits[f"seg{segMatch}"] ## check if segment is marked as present
                    fc = 'k' if segInvolved else 'w' ## black if segment traveled along edge, white otherwise
                    
                    ax.add_patch(Rectangle((sx, y - segHeight / 2), width=segWidth, height=segHeight, fc=fc, ec='k', clip_on=False))
                    
                    if i == (len(segOrder) - 1): ## last (horizontally) rectangle, label with number
                        ax.text(sx + segWidth * 1.05, y, f"#{reassortmentCounter}", ha='left', va='center', size=12)
                    if reassortmentCounter == 1: ## first (top-most) reassortment, label segment
                        ax.text(sx + segWidth / 2, y + segHeight, seg, ha='left', va='center', rotation=90, rotation_mode='anchor', size=12)

    ax.add_collection(LineCollection(pointerLine, ls='--', lw=1, color=pointerLineColours, clip_on=False)) ## pointer line that connects end of reticulation to schematic of segment involvement matrix
    ax.add_collection(LineCollection(branches, ls='--', lw=2, color=colours)) ## plot reticulation (originates as normal branch, heads down to destination branch
    
    return ax

def plot_tree_matrix(treeAx, matrixAx, tree, labelDict, colourDict=None, columnOrder=None, width=1.0, **kwargs):

    from matplotlib.patches import Rectangle
    from matplotlib.collections import PatchCollection
    
    tree.plot_tree(treeAx)

    if columnOrder is None:
        columnOrder = sorted(labelDict.keys())
    else:
        assert len(labelDict) == len(columnOrder), f"labelDict has length {len(labelDict)} but columnOrder has length {len(columnOrder)}"
    localKwargs = dict(kwargs)
    
    if colourDict is not None: ## colourDict was provided
        for col in columnOrder: ## check that column colours are provided
            assert col in colourDict, f"Column {col} not found in colourDict"
    
    cells = []
    colours = []
    for i, col in enumerate(columnOrder): ## iterate over matrix columns
        x = (i * width)

        for j, taxon in enumerate(tree.get_external()): ## iterate over taxa in tree
            y = taxon.y - 0.5

            if taxon.name in labelDict[col]: ## have label for taxon
                label = labelDict[col][taxon.name]
            else:
                label = '' ## blank label, warn user
                logger.warning(f"labelDict missing taxon {taxon.name} at column {col}")

            if label in colourDict[col]: ## colour matrix if label available
                fc = colourDict[col][label]
            else:
                fc = 'none'
                logger.warning(f"Taxon {taxon.name} label '{label}' not assigned a colour, defaulting to invisible")
            
            cell = Rectangle((x, y), width=width, height=1.0)
            cells.append(cell)
            colours.append(fc)
    
    assert 'fc' not in localKwargs and 'facecolor' not in localKwargs, f"Not allowed to set facecolor of matrix cells since they're set via colourDict."
    matrixAx.add_collection(PatchCollection(cells, color=colours, **localKwargs))

    treeAx.set_ylim(0, tree.ySpan + 0.5)
    matrixAx.set_ylim(treeAx.get_ylim())

    matrixStart = 0
    matrixEnd = (len(columnOrder) * width)
    xticks = np.linspace(matrixStart + width/2, matrixEnd - width/2, len(columnOrder))
    matrixAx.set_xticks(xticks)
    matrixAx.set_xticklabels(columnOrder, rotation=90)
    matrixAx.set_yticks([])
    matrixAx.set_yticklabels([])
    
    matrixAx.set_xlim(matrixStart, matrixEnd)

    return treeAx, matrixAx
##############

def _compute_consensus(alnFile, SNPs=None, validNucleotideFxn=None, alnFmt='fasta'):
    """
    Computes consensus from an alignment file provided.
    Optionally can compute SNP consensus
    """
    from collections import Counter
    from Bio import SeqIO
    
    alnDict = SeqIO.to_dict(SeqIO.parse(alnFile, alnFmt))
    
    if validNucleotideFxn is None: validNucleotideFxn = lambda nt: nt.upper() in ['A','C','T','U','G','-']

    if SNPs is None:
        positions = range(len(alnDict[list(alnDict.keys())[0]])) ## get length of first sequence == alignment length
    else:
        assert isinstance(SNPs, list), f"Parameter SNPs is of type {type(SNPs)}, not list"
        positions = SNPs

    consensus = []
    for i in positions: ## iterate over positions (either range or list of pre-computed variable sites)

        clean_column = [alnDict[seq][i] for seq in alnDict if validNucleotideFxn(alnDict[seq][i].upper())] ## filter to valid nucleotides
        assert len(clean_column) > 0, f"Position {i+1} (1-based indexing) contains no valid nucleotides"
        consensusNt = clean_column[0] if len(set(clean_column)) == 1 else Counter(clean_column).most_common(1)[0][0]
        
        consensus.append(consensusNt)
    
    return ''.join(consensus)

def _get_refSeq(refSeq, validNucleotideFxn=None, alnFile=None, alnFmt=None, refSeqFmt=None):
    """
    Retrieve reference sequence - compute consensus from alignment, fetch from alignment file or grab from path.
    Returns sequence as str
    """
    from Bio import SeqIO
    
    summarySeq = None
    if refSeq == 'consensus': ## compute consensus
        assert alnFile is not None and validNucleotideFxn is not None and alnFmt is not None, f"Computing consensus requires alnFile, alnFmt, validNucleotidesFxn to be set."
        if refSeqFmt:
            logger.warning(f"Using consensus sequence as reference, refSeqFmt parameter ({refSeqFmt}) will be ignored.")
        summarySeq = _compute_consensus(alnFile=alnFile, validNucleotideFxn=validNucleotideFxn, alnFmt=alnFmt)
    elif os.path.exists(refSeq):
        assert refSeqFmt is not None, f"Loading an external reference requires refSeqFmt to be set."
        summarySeq = str(SeqIO.read(refSeq, refSeqFmt).seq)
    else:
        assert alnFmt is not None, f"Identifying a reference in alignment requires alnFmt to be set."
        alnDict = SeqIO.to_dict(SeqIO.parse(alnFile, alnFmt)) ## load sequence
        assert refSeq in alnDict, f"refSeq {refSeq} not found in alignment {alnFile} provided."
        if refSeqFmt:
            logger.warning(f"Using sequence in alnDict as reference, refSeqFmt parameter ({refSeqFmt}) will be ignored.")
        summarySeq = str(alnDict[refSeq].seq) ## grab reference sequence without filtering invalid nucleotides

    if summarySeq is None:
        raise Exception(f"Could not retrieve reference sequence.")
    
    return summarySeq

def _default_variable_site_selection_fxn(columnDict, validNtFxn, refNt):
    """
    Takes a dict of {seq name: nucleotide}, validNtFxn to identify whether a nucleotide is valid, and nucleotide at reference sequence
    This function identifies whether the column dict contains variable sites (at least 2 unique values) and that the 2nd most common allele is found it at least 2 sequences.
    """
    variable = False
    clean_column = [columnDict[seq] for seq in columnDict if validNtFxn(columnDict[seq])] ## filter to valid nucleotides
    
    if len(set(clean_column)) >= 2 and Counter(clean_column).most_common()[1][1] >= 2: ## column must be variable AND second most common element is found more than once
        variable = True

    return variable

def get_variable_aln_sites(alnFile, refSeq='consensus', selectionFxn=None, validNucleotideFxn=None, alnFmt='fasta', refSeqFmt=None, trimStart=0, trimEnd=None):
    """
    Takes an alignment file path.
    refSeq is consensus (computed from alignment provided), a specific sequence in the alignment or an external file.
    selectionFxn is a function that takes: a dict of {seq name: nucleotide} at a given alignment column, default: _default_variable_site_selection_fxn that filters column down to valid nucleotides, and nucleotide at reference sequence; default: True for columns where the 2nd most common nt is found in >=2 seqs
    validNucleotideFxn is a function that determines whether a nucleotide is valid (default is nt == A, C, T, U, G or -)
    alnFmt is sequence format of alnFile (default: 'fasta')
    refSeqFmt is sequence format of external reference file (default: 'fasta')
    trimStart, trimEnd is how many positions to ignore at either end of the alnFile
    returns indices of sites that pass selectionFxn
    """
    from collections import Counter
    from Bio import SeqIO
    
    alnDict = SeqIO.to_dict(SeqIO.parse(alnFile, alnFmt)) ## load sequence
    
    if validNucleotideFxn is None: validNucleotideFxn = lambda nt: nt.upper() in ['A','C','T','U','G','-']
    
    if selectionFxn is None: selectionFxn = _default_variable_site_selection_fxn
    
    if trimStart is None: trimStart = 0

    alnLs = [len(seq) for seq in alnDict.values()]
    assert len(set(alnLs)) == 1, f"Sequences in alignment have different lengths implying they're not aligned. Unique sequence lengths in alignment: {set(alnLs)}."
    alnL = alnLs[0] ## all seqs have same length, get first seq length
    
    if trimEnd is None:
        trimEnd = alnL
    else:
        trimEnd = alnL - trimEnd

    assert refSeq == 'consensus' or refSeq in alnDict or os.path.exists(refSeq), f"refSeq {refSeq} is not 'consensus', a sequence available in alnDict or a recognised path."
    if os.path.exists(refSeq) and refSeqFmt is None: ## external reference without provided format, assume fasta
        refSeqFmt = 'fasta'
    summarySeq = _get_refSeq(refSeq=refSeq, alnFile=alnFile, alnFmt=alnFmt, refSeqFmt=refSeqFmt, validNucleotideFxn=validNucleotideFxn)
    
    assert len(summarySeq) == alnL, f"External reference sequence length {len(summarySeq)} does not match alignment sequence length"
    
    columns = [{s: str(alnDict[s].seq)[i] for s in alnDict} for i in range(alnL)] ## extract every column as a {seq: nt} dict
    SNPs = [i for i in range(alnL) if selectionFxn(columns[i], validNucleotideFxn, summarySeq[i]) and trimStart <= i <= trimEnd] ## identify indices that satisfy filterFxn and are within desired range
    
    return SNPs

def plot_snp_alignment(alnAx, SNPs, alnFile, tree, refSeq='consensus', ntColours=None, textKwargs={}, 
                       rectangleKwargs={}, treeAx=None, fmtSeqNamesFxn=None, treeKwargs={}, alnFmt='fasta', 
                       validNucleotideFxn=None, coding=False, gffFile=None, featType=None, geneName=None, refSeqFmt=None, 
                       plotORFs=False, minFeatLen=1000, offsetORFs=0.1):
    """
    takes alnAx axes object where the alignment will be plotted
    SNPs is a list of variable sites to be included (0-indexed)
    alnFile is the alignment file
    tree is the tree, used to sort alignment in tree y-axis order; ## NOTE - should assert that number of sequences matches number of tips in tree
    refSeq is the reference sequence to be used to highlight sequence differences:
    - if 'consensus', then a consensus sequence will be computed from the alignment provided and used as a reference
    - if refSeq is a sequence in the alignment, that sequence becomes the reference
    - if refSeq is a path, it's assumed that an external reference file was provided (must contain a single record!)
    ntColours is a dict that maps every nucleotide state to a colour; note, don't forget ambiguous state (N, K, M, R, etc.) colours if using your own
    textKwargs is a kwargs dict for text labels at changed nucleotides (or reference)
    rectangleKwargs is a kwargs dict for rectangles that represent each nucleotide in each sequence
    treeAx is an axes object for placing a tree; can be None if you don't want the tree plotted
    fmtSeqNamesFxn is a function that takes sequence name from a tree and returns a modified version (e.g. you want to remove extraneous information)
    treeKwargs is a kwargs dict used when plotting a tree (e.g. you want to change branch colours)
    alnFmt is the format of the alignment file (default: 'fasta')
    validNucleotideFxn is a function that filters alignment columns down to valid states; default {A, C, T, U, G, -} are considered as valid 
    coding is a Boolean for whether the alignment is coding or not; if a GFF file is provided the features in the file will be used to infer aa changes, otherwise it's assumed the entire alignment is a single coding sequence
    gffFile is a path (or handle) to a GFF file that defines sequence features; ## NOTE - should assert that alignment length equals 'nuc' feature
    featType is a string used for fetching features from the GFF file; if not specified assumed to be 'gene'
    geneName is a string used for getting feature names from the GFF file; if not specified assumed to be 'gene_name'
    refSeqFmt is the format of an external reference file (i.e. if refSeq is a path)
    plotORFs is a Boolean for whether to plot a schematic diagram of the alignment with feature annotations; errors if set to True without a GFF file
    minFeatLen is minimum feature length for adding text to sequence schematic
    offsetORFs is a float that represents a fraction of the tree's y-axis height to be used for positioning GFF file features below the alignment
    """ 
    from matplotlib.patches import Rectangle
    from matplotlib.collections import PatchCollection, LineCollection
    from Bio import SeqIO
    
    alnDict = SeqIO.to_dict(SeqIO.parse(alnFile, alnFmt))
    
    height = 0.95
    width = 1.0
    localTextKwargs = dict(textKwargs)
    localRectangleKwargs = dict(rectangleKwargs)
    localTreeKwargs = dict(treeKwargs)
    
    if 'color' not in localTextKwargs: localTextKwargs['color'] = 'k'
    if 'size' not in localTextKwargs: localTextKwargs['size'] = 10
    if 'zorder' not in localTextKwargs: localTextKwargs['zorder'] = 100
    if 'ha' not in localTextKwargs and 'horizontalalignment' not in localTextKwargs: localTextKwargs['ha'] = 'center'
    if 'va' not in localTextKwargs and 'verticalalignment' not in localTextKwargs: localTextKwargs['va'] = 'center'
    
    if ntColours is None: ## default nt colours
        ntColours={'A': '#D0694A', 'C': '#77BEDB', 'T': '#48A365', 'U': '#48A365', 'G': '#E1C72F', 
                 '-': 'white', 'N':'dimgrey', 
                 'K': 'dimgrey', 'Y': 'dimgrey', 'M': 'dimgrey', 'W': 'dimgrey', 'R': 'dimgrey'}
    
    ###### check refSeq is valid, grab reference sequence accordingly
    assert refSeq == 'consensus' or refSeq in alnDict or os.path.exists(refSeq), f"refSeq {refSeq} is not 'consensus', a sequence available in alnDict or a recognised path."
    if validNucleotideFxn is None: validNucleotideFxn = lambda nt: nt in ['A','C','T','U','G','-']
    summarySeq = _get_refSeq(refSeq=refSeq, alnFile=alnFile, alnFmt=alnFmt, refSeqFmt=refSeqFmt, validNucleotideFxn=validNucleotideFxn)
    ######
    if gffFile:
        seqFeatures = _load_gff(gffFile)
        if featType is None: featType = 'gene'
        if geneName is None: geneName = 'gene_name'
    else:
        seqFeatures = None
        if coding:
            logger.warning("Alignment is coding but no GFF file provided.")
    ###### compute xtick positions for alignment columns
    window = 3
    storeSite = SNPs[0]
    xticks = []
    cumulativeX = -1
    for i, pos in enumerate(SNPs):
        if storeSite + window < pos: ## next site is beyond window
            cumulativeX += 1 + (pos - storeSite) * 0.0002 ## add a bit of extra space once we're far enough away
        else:
            cumulativeX += 1

        xticks.append(cumulativeX)
        storeSite = pos
    ###########
    patches = []
    patchFaceColours = []
    patchEdgeColours = []
    
    for k in sorted(tree.get_external(), key = lambda leaf: leaf.y):
        assert k.name in alnDict, f"Tip name {k.name} not present in alnDict."
        
        y = k.y - 0.5

        for i, pos in enumerate(SNPs):
            x = xticks[i]
            
            curNt = str(alnDict[k.name][pos]).upper() ## grab current sequence's nt
            refNt = summarySeq[pos] ## grab consensus (at SNP index) or reference sequence nucleotide (at entire seq index)
            
            nt_block = Rectangle((x, y), width = width, height = height)
            
            if refNt != curNt: ## not matching with reference sequence - highlight, add text
                fc = ntColours[curNt]
                ec = 'w'
                alnAx.text(x + 0.5, y + 0.5, curNt, **localTextKwargs)
            elif refNt == '-': ## reference sequence has a gap - there's insertions elsewhere, so no marker
                fc = 'none' ## invisible block
                ec = 'none'
            else: ## match with reference/consensus
                fc = 'lightgray'
                ec = 'none'
                if k.name == refSeq: ## add text for reference seq
                    alnAx.text(x + 0.5, y + 0.5, curNt, **localTextKwargs)
            
            patches.append(nt_block)
            patchFaceColours.append(fc)
            patchEdgeColours.append(ec)

        if refSeq in alnDict and refSeq == k.name:
            alnAx.add_patch(Rectangle((0, k.y - height/2), width=max(xticks) + width, height=height, facecolor='none', edgecolor='k', lw=1, zorder=10000, clip_on=False)) ## outline reference
    
    ###### plot consensus or external reference sequence if provided underneath alignment
    if refSeq == 'consensus' or os.path.exists(refSeq):
        y = -1
        for i, pos in enumerate(SNPs):
            x = i
            nt_block = Rectangle((xticks[i], y), width = width, height = height)
            
            alnAx.text(xticks[i] + 0.5, y + 0.5, summarySeq[pos].upper(), color = 'k', size = 10, ha = 'center', va = 'center')
            
            patches.append(nt_block)
            patchFaceColours.append('lightgray')
            patchEdgeColours.append(ec)

        alnAx.add_patch(Rectangle((0, y), width=max(xticks) + width, height=height, facecolor='none', edgecolor='k', lw=1, zorder=10000, clip_on=False)) ## outline reference
    ############ dump all rectangles, adjust axes
    alnAx.add_collection(PatchCollection(patches, facecolors = patchFaceColours, edgecolors = patchEdgeColours))

    ######
    alnAx.xaxis.tick_top()
    alnAx.set_xticks([x + width/2 for x in xticks])

    xTickLabels = []

    for i, pos in enumerate(SNPs):
        if coding:
            fmtPosition = _format_coding_aln_column(site=pos, alnDict=alnDict, referenceSeq=summarySeq, seqFeatures=seqFeatures, featType=featType)
        else:
            fmtPosition = pos + 1 ## convert to 1-indexed

        xTickLabels.append(fmtPosition)
    
    if len(xTickLabels) == 0: print(f"xTickLabels not set")
    alnAx.set_xticklabels(xTickLabels, rotation = 90)
    alnAx.tick_params(size = 0)

    alnBottom = -1 if (refSeq == 'consensus' or os.path.exists(refSeq)) else 0
    ######## plotting representative sequence features
    if plotORFs:

        # adjust x-axis shoulder position
        repulsionStrength = 0.01
        xs = np.array(xticks, float)
        for _ in range(12):  # 10 iterations of simple repulsion
            diffs = xs[:,None] - xs[None,:]
            xs += repulsionStrength * np.tanh(1 / (diffs + 1e-6)).sum(axis=1)
        
        spread = dict(zip(SNPs, xs))
        ###
        
        assert gffFile, f"Must provide gffFile for plotting ORFs."
        ORFwidth = min([tree.ySpan * 0.05, 1])
        ORFwidth = max([ORFwidth, tree.ySpan * 0.05])
        rescaleORFs = len(summarySeq) / (max(xticks) + width)
        offsetOrfY = tree.ySpan * offsetORFs ## y-axis coordinate where genome sits, offset by a fraction of tree length
        maxORFheight = plot_seq_features(alnAx, gffFile=gffFile, xy=(0, -offsetOrfY), width=ORFwidth, rescale=rescaleORFs, geneName='gene_name', minFeatLen=minFeatLen)
        alnAx.eventplot([i / rescaleORFs for i in SNPs], lineoffsets=-offsetOrfY, linelengths=ORFwidth*0.4, color='k', clip_on=False)
        yBottom = alnBottom - offsetOrfY - maxORFheight - ORFwidth * 0.5 ## adjust bottom to show ORFs starting for bottom of alignment
        
        connections = []
        for i, pos in enumerate(SNPs):
            shoulder_x = spread[pos] #/ rescaleORFs
            connections.append(((xticks[i] + width / 2, alnBottom), ## bottom middle of alignment column
                                (shoulder_x + width / 2, alnBottom - offsetOrfY * 0.3 + ORFwidth * 0.2), ## retain x position, move halfway down y-axis
                                (pos / rescaleORFs, -offsetOrfY + ORFwidth * 0.2))) ## connect to sequence position
        
        alnAx.add_collection(LineCollection(connections, color='gray', lw=0.5, clip_on=False))
        if coding == False:
            logger.warning(f"plotORFs == True, but coding == False; Alignment columns will not be formatted with amino acid changes to ORFs.")
    else: ## not plotting ORFs, bottom of plot is where alignment ends
        yBottom = alnBottom
    ###### tree part
    if treeAx:
        tree.plot_tree(treeAx, autoSort=False, **localTreeKwargs)
        
        for k in tree.get_external():
            treeAx.plot([k.x, tree.treeHeight * 1.01], [k.y, k.y], color = 'gray', ls = '--', lw = 0.5)

        treeAx.xaxis.tick_top()
        treeAx.tick_params(rotation=90)
        treeAx.set_xlim(-tree.treeHeight * 0.01, max(tree.get_parameter_list('x')) * 1.01)
        treeAx.set_ylim(yBottom - 0.05, tree.ySpan + 0.05)

        from baltic.bt_utils import clean_axes
        clean_axes(treeAx, hideSpines = ['bottom', 'left', 'right'], removeTickLabels='y')
    #######
    if fmtSeqNamesFxn is None: fmtSeqNamesFxn = lambda k: f"{k.name}"
    yticks = [k.y for k in tree.get_external()]
    yTickLabels = [fmtSeqNamesFxn(k) for k in tree.get_external()]

    if refSeq == 'consensus':
        yticks.insert(0, -height/2)
        yTickLabels.insert(0, 'consensus sequence')
    elif os.path.exists(refSeq):
        yticks.insert(0, -height/2)
        yTickLabels.insert(0, 'reference sequence')
    
    alnAx.yaxis.tick_right()
    alnAx.set_yticks(yticks)
    alnAx.set_yticklabels(yTickLabels)
    #######
    
    alnAx.set_ylim(yBottom - 0.05, tree.ySpan + 0.05)
    alnAx.set_xlim(-0.1, max(xticks) + width + 0.05)
    [alnAx.spines[loc].set_visible(False) for loc in alnAx.spines]
    
    return alnAx

def _identify_gene(site, seqFeatures, featType, geneName):
    """
    Identifies which feature (if any) a given sequence position is in and the index of where the feature begins.
    If sequence position is within a feature, returns feature name and its starting position; otherwise returns None, None 
    Site is 0-indexed position in alignment.
    seqFeatures is a list of Biopython seq features (assumed loaded from a GFF file earlier)
    featType and geneName are the feature type and feature key for extraction from GFF file
    """
    hits = []

    for feat in seqFeatures:
        start = int(feat.location.start)
        end   = int(feat.location.end)

        if start <= site < end and feat.type == featType:
            name = feat.qualifiers.get(geneName, ['?'])[0]
            geneStart = site - start
            hits.append((name, geneStart))

    return hits   # may be empty, 1, or many


def _format_coding_aln_column(site, alnDict, referenceSeq, seqFeatures=None, featType='gene', geneName='gene_name'):
    """
    site: 0-indexed alignment column
    alnDict: {name: SeqRecord} dict
    referenceSeq: str (same length as alignment)
    seqFeatures: optional list of SeqFeatures from GFF
    featType and geneName are the feature type (default 'gene') and feature key (default 'gene_name') for extraction from GFF file
    returns the label for a given alignment column, e.g.:
    - "1410 nt" if site was outside of a coding feature
    - "23012 nt S: E484K (A) / E484A (C)" for a site that changes the amino acid sequence (lists all aa variants and nucleotide change associated with each outcome)
    - "24853 nt S: 1097 aa" for a site that's within a coding feature but the mutation is synonymous
    """
    from Bio.Seq import Seq

    # Identify all overlapping ORFs
    hits = []
    if seqFeatures is not None:
        hits = _identify_gene(site, seqFeatures, featType, geneName)

    # Outside coding — no ORF overlaps
    if not hits:
        return f"{site + 1} nt"

    annotations = []

    for gene, geneStart in hits:

        # Compute codon position within this ORF
        aaSite   = geneStart // 3
        codonPos = geneStart % 3
        codonStart = site - codonPos

        refCodon = referenceSeq[codonStart: codonStart + 3].upper()

        if len(refCodon) != 3:
            annotations.append(f"{gene}: {site+1} nt")
            continue

        aaRef = "-" if refCodon == "---" else str(Seq(refCodon).translate())

        # Collect codons from all sequences
        allCodons = [
            str(alnDict[name].seq[codonStart: codonStart + 3]).upper()
            for name in alnDict
        ]
        varCodons = {c for c in allCodons if len(c) == 3 and c[codonPos] != refCodon[codonPos]}

        fmtAAsite = aaSite + 1
        aaChanges = {}

        for altCodon in varCodons:
            altNt = altCodon[codonPos]

            if 0 < altCodon.count("-") < 3: ## 1 or 2 gaps in codon - frameshift
                aaNew = "frmShft"
            elif altCodon == "---":
                aaNew = "-"
            else:
                aaNew = str(Seq(altCodon).translate())

            if aaNew == "X":
                continue

            fmtRefAA = r"$\Delta$" if aaRef == "-" else aaRef
            fmtNewAA = r"$\Delta$" if aaNew == "-" else aaNew

            if fmtRefAA == fmtNewAA:
                continue

            key = (fmtRefAA, fmtAAsite, fmtNewAA)
            aaChanges.setdefault(key, set()).add(altNt)

        if aaChanges:
            effects = " / ".join(
                f"{ref}{pos}{new} ({','.join(sorted(ntSet))})"
                for (ref, pos, new), ntSet in sorted(aaChanges.items())
            )
        else:
            effects = f"{fmtAAsite} aa"  # synonymous

        annotations.append(f"{gene}: {effects}")

    # Combine ORF annotations, separated by " | "
    return f"{site + 1} nt " + " | ".join(annotations)


def _load_gff(gffFile):
    """
    Loads and returns features from a GFF file.
    """
    try:
        from BCBio import GFF
    except ImportError:
        raise ImportError(f"BCBio-GFF needs to be installed for processing GFF files. Install with: pip install bcbio-gff")

    handle = open(gffFile, 'r') if isinstance(gffFile, str) else gffFile

    features = []
    for rec in GFF.parse(gffFile):
        features.extend(rec.features)

    return features

def _assign_tracks(features):
    """
    Assign minimal y-level tracks to features such that overlapping CDSs 
    do not share a track.

    features: list of (feat_obj, start, end)
    Returns: dict {feature_index: track_index}
    """
    # Give each feature a stable index
    indexed = [(i, feat, start, end) for i, (feat, start, end) in enumerate(features)]
    # Sort by start coordinate
    indexed.sort(key=lambda x: x[2])

    tracks = []   # list of lists of (start,end)
    assignment = {}

    for idx, feat, start, end in indexed:
        placed = False

        # Try existing tracks
        for t_idx, intervals in enumerate(tracks):
            if all(end <= s or start >= e for (s, e) in intervals):
                intervals.append((start, end))
                assignment[idx] = t_idx
                placed = True
                break

        if not placed:
            tracks.append([(start, end)])
            assignment[idx] = len(tracks) - 1

    return assignment

def plot_seq_features(ax, gffFile, xy=None, rescale=None, width=None, geneName='gene_name', arrowKwargsFxn=None, textArgsFxn=None, textKwargsFxn=None, minFeatLen=1000):
    """
    Takes axes and GFF file and plots arrows
    xy is position where annotations start
    rescale is a parameter that controls the length of annotation - if aln is 100nt long, rescale set to 100.0 will plot all annotations within interval [0.0, 1.0]
    geneName is the field name of a feature in the GFF file
    arrowKwargsFxn is a function that accepts the name of a feature (e.g. "ORF1a", "S", etc.) and creates a kwargs dict to be used when plotting feature arrows
    textKwargsFxn is a function that accepts the name of a feature (e.g. "ORF1a", "S", etc.) and returns a kwargs dict to be used to modify text that's plotted
    textArgsFxn is a function that accepts the name of a feature (e.g. "ORF1a", "S", etc.) and returns an args dict to be used for adding text
    minFeatLen is a threshold for adding text to features - no text is plotted for features less than this length
    """
    from matplotlib.patches import FancyArrow
    
    if xy is None: xy = (0, 0)
    if rescale is None: rescale = 1.0
    if width is None: width = 1.0
    
    start_x, start_y = xy
    
    features = _load_gff(gffFile) ## load GFF annotations
    
    cds_intervals = []
    for feat in features:
        start = int(feat.location.start)
        end = int(feat.location.end)
        cds_intervals.append((feat, start, end))
    
    track_assignment = _assign_tracks(cds_intervals) ## position ORFs so they don't overlap
    
    ys = []
    for i, (feat, start, end) in enumerate(cds_intervals):

        rescaled_start = start / rescale
        rescaled_end = end / rescale
        
        name = feat.qualifiers.get(geneName, ['?'])[0]
        track = -track_assignment[i]

        length = end - start
        rescaled_length = length / rescale
        
        head_length = 300 if length > 300 else max(20, int(length * 0.3))
        rescaled_head_length = head_length / rescale
        
        y = start_y - width * 0.1 + (track * width)
        ys.append(y)
        if name == 'nuc':
            # Genome backbone on track 0
            ax.plot([start_x + rescaled_start, start_x + rescaled_end], [start_y, start_y], color='k', lw=3, zorder=0, clip_on=False)
            ax.text(start_x + rescaled_end * 1.01, start_y, f"{end}", ha='left', va='center', color='k', clip_on=False)
            continue
    
        ##############
        defaultArrowKwargs = {'width': width, 'head_width': width, 'head_length': rescaled_head_length, 
                              'length_includes_head': True, 'facecolor': 'lightgray', 'edgecolor': 'k', 
                              'zorder': 10}
        
        localArrowKwargs = dict(defaultArrowKwargs) ## copy defaults
        
        if arrowKwargsFxn is not None:
            arrowFxnDict = arrowKwargsFxn(name) ## grab dict returned by function
            if 'fc' in arrowFxnDict: arrowFxnDict['facecolor'] = arrowFxnDict['fc']; arrowFxnDict.pop('fc')
            if 'ec' in arrowFxnDict: arrowFxnDict['edgecolor'] = arrowFxnDict['ec']; arrowFxnDict.pop('ec')
            
            for key in arrowFxnDict: ## overwrite values
                localArrowKwargs[key] = arrowFxnDict[key]
        
        arrow = FancyArrow(x = start_x + rescaled_start, y = y, 
                           dx = rescaled_length, dy = 0, **localArrowKwargs)
        ax.add_patch(arrow)
        ###############
        defaultTextArgs = {'x': start_x + rescaled_start + rescaled_length/2, 'y': y , 's': name}
        localTextArgs = dict(defaultTextArgs)
        
        if textArgsFxn is not None:
            textFxnDict = textArgsFxn(name) ## grab dict returned by function
            if 'fc' in textFxnDict: textFxnDict['facecolor'] = textFxnDict['fc']; textFxnDict.pop('fc')
            if 'ec' in textFxnDict: textFxnDict['edgecolor'] = textFxnDict['ec']; textFxnDict.pop('ec')

            for key in textFxnDict: ## overwrite values
                localTextArgs[key] = textFxnDict[key]
        #############
        defaultTextKwargs = {'ha': 'center', 'va': 'center', 'zorder': 100, 'rotation': 0, 'rotation_mode': 'anchor'}
        localTextKwargs = dict(defaultTextKwargs)
        
        if textKwargsFxn is not None:
            textFxnDict = textKwargsFxn(name) ## grab dict returned by function
            if 'fc' in textFxnDict: textFxnDict['facecolor'] = textFxnDict['fc']; textFxnDict.pop('fc')
            if 'ec' in textFxnDict: textFxnDict['edgecolor'] = textFxnDict['ec']; textFxnDict.pop('ec')

            for key in textFxnDict: ## overwrite values
                localTextKwargs[key] = textFxnDict[key]

        if minFeatLen <= length:
            ax.text(**localTextArgs, **localTextKwargs)

    maxHeight = abs(min(ys) - start_y)
    
    return maxHeight