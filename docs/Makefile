# Minimal makefile for building the docs with Sphinx
SPHINXBUILD   = sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

.PHONY: html clean

html:
	$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" -q

clean:
	rm -rf $(BUILDDIR)
