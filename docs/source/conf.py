# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'baltic'
copyright = '2026, Gytis Dudas & Barney Potter'
author = 'Gytis Dudas & Barney Potter'
release = 'v0.9 (Birch)'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'matplotlib': ('https://matplotlib.org/stable/', None),
    'biopython': ('https://biopython.org/docs/latest/', None),
}

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'sphinx_rtd_theme'
# html_theme = 'furo'
html_theme = 'alabaster'
html_static_path = ['_static']
html_css_files = ['copybutton.css']
html_js_files = ['copybutton.js']

# Let ``furo`` use its built-in sidebar configuration.
