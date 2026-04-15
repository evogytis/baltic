# Minimal makefile for Sphinx documentation
#

# Prefer the repo-managed docs environment when it exists so local builds use
# the same theme and dependency set declared in anaconda-project.yml.
LOCALSPHINX   := $(CURDIR)/envs/default/bin/sphinx-build
SPHINXOPTS    ?=
ifeq ($(wildcard $(LOCALSPHINX)),)
SPHINXBUILD   ?= sphinx-build
else
SPHINXBUILD   ?= $(LOCALSPHINX)
endif
SOURCEDIR     = source
BUILDDIR      = build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
