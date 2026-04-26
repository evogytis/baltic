[![Build Status](https://travis-ci.com/evogytis/baltic.svg?branch=master)](https://travis-ci.com/evogytis/baltic)
[![downloads](https://anaconda.org/bioconda/baltic/badges/downloads.svg)](https://anaconda.org/bioconda/baltic)

# baltic

`baltic` is a lightweight Python library for parsing, manipulating, and visualizing phylogenetic trees. It supports Newick, Nexus, and Nextstrain/Auspice JSON inputs and provides plotting helpers built on top of `matplotlib`.

The name expands to Backronymed Adaptable Lightweight Tree Import Code.

## Installation

Install from PyPI with `pip`:

```bash
pip install baltic
```

`baltic` depends on `numpy` and `matplotlib`. If you plan to load JSON trees from URLs, install `requests` as well.

For local development and documentation builds, create the single conda environment declared in [`baltic.yaml`](/Users/bip7/Tools/baltic/baltic.yaml), then install the package into it:

```bash
conda env create -f baltic.yaml
conda activate baltic
pip install -e . --no-build-isolation
```

## Quick start

### Parse a tree string

```python
import baltic as bt

tree_string = "((A:1.0,B:2.0):1.0,C:3.0);"
ll = bt.make_tree(tree_string, treeType="divergence")

ll.treeStats()
```

### Load trees from files

The public loader functions use `snake_case` names:

```python
import baltic as bt

newick_tree = bt.load_newick("tree.nwk", treeType="divergence")
nexus_tree = bt.load_nexus("tree.nex", treeType="time")
json_tree, json_meta = bt.load_JSON(
    "https://nextstrain.org/charon/getDataset?prefix=/dengue/denv1",
    treeType="time",
)
```

Common options:

- `treeType="divergence"` for branch lengths in relative units.
- `treeType="time"` for time-calibrated trees.
- `absoluteTime=True` when tip dates should be parsed from tip labels.
- `tipRegex` and `dateFmt` to control how sampling dates are extracted.

For example, if tip names end with ISO dates:

```python
time_tree = bt.load_nexus(
    "example.tree",
    treeType="time",
    tipRegex=r"\|([0-9\-]+)$",
    dateFmt="%Y-%m-%d",
    absoluteTime=True,
)
```

## Basic tree operations

```python
tips = ll.get_external()
internal = ll.get_internal()

stats = ll.treeStatsDict()
subtree = ll.subtree(tips[0])
```

Useful methods on `Tree` objects include:

- `traverse_tree()` to recompute heights and cached traversal state.
- `sort_branches()` to order branches for plotting.
- `rename_tips()` to apply a translation table.
- `set_absolute_time()` to assign decimal dates on time trees.
- `plot_tree()`, `plot_points()`, and `plot_text()` for visualization.

## Plotting example

```python
import matplotlib.pyplot as plt
import baltic as bt

ll = bt.load_newick("tree.nwk", treeType="divergence")

fig, ax = plt.subplots(figsize=(8, 10))
ll.plot_tree(ax)
ll.plot_text(ax, targetFxn=lambda k: k.is_leaf())

ax.set_axis_off()
fig.tight_layout()
plt.show()
```

`plot_tree()` supports rectangular, circular, and unrooted layouts through the `treeType` plotting argument.

## Documentation

The Sphinx sources live in [`docs/source`](/Users/bip7/Tools/baltic/docs/source).

Build the docs from the same `baltic` conda environment:

```bash
cd docs
make html
```

## Development

Run the test suite with:

```bash
python tests/testsuite.py
```

## License

Copyright 2016 [Gytis Dudas](https://twitter.com/evogytis).

Licensed under the [GNU GPL v3.0](https://github.com/evogytis/baltic/blob/master/LICENSE).
