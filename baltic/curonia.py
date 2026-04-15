import warnings
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize
from baltic.bt_utils import calendar_to_decimal_date


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
