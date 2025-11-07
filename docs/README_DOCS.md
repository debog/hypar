# HyPar Documentation

This directory contains the Read the Docs documentation for HyPar.

## Documentation Structure

- **index.md** - Main landing page with project overview
- **installation.md** - Detailed build and installation instructions
- **usage.md** - Complete guide to input files and running simulations
- **numerical_methods.md** - Description of numerical algorithms
- **developer_guide.md** - Guide for adding new physical models
- **api.md** - API reference (stub for autodoc)
- **conf.py** - Sphinx configuration

## Building Documentation Locally

### Prerequisites

Install the required Python packages:

```bash
pip install -r requirements-docs.txt
```

This installs:
- Sphinx (documentation generator)
- MyST Parser (Markdown support)
- sphinx-rtd-theme (Read the Docs theme)
- sphinx-autodoc-typehints (for API docs)

### Build HTML Documentation

From the `docs/` directory:

```bash
make html
```

Then open `_build/html/index.html` in your web browser.

### Clean Build

To start fresh:

```bash
make clean
make html
```

## Documentation Content

### For Users

1. **Installation Guide** - Comprehensive instructions for:
   - Basic serial/MPI builds
   - Optional dependencies (PETSc, CUDA, FFTW, etc.)
   - Troubleshooting common issues

2. **Usage Guide** - Complete reference for:
   - All input file formats (solver.inp, boundary.inp, physics.inp, etc.)
   - Running simulations (serial and parallel)
   - Output files and interpretation
   - Quick start example

3. **Numerical Methods** - Detailed explanation of:
   - Spatial discretization schemes (WENO, MUSCL, etc.)
   - Time integration methods (explicit RK, implicit/IMEX via PETSc)
   - Upwinding schemes
   - Stability and CFL conditions

### For Developers

1. **Developer Guide** - Instructions for:
   - Code structure overview
   - Adding new physical models (step-by-step)
   - Required and optional functions
   - Build system integration
   - Testing and validation

2. **API Reference** - Code documentation (requires autodoc setup)

## Content Sources

Documentation was compiled from:
- `doc/Input_Files.md` → `usage.md`
- `doc/Numerical_Method.md` → `numerical_methods.md`
- `doc/Adding_Physical_Models.md` → `developer_guide.md`
- `configure.ac` + README → `installation.md`
- New comprehensive overview → `index.md`

## Publishing to Read the Docs

This documentation is designed to be automatically built by Read the Docs when:
1. The `.readthedocs.yaml` configuration file is present (already created)
2. The repository is linked to a Read the Docs project

Configuration:
- Python 3.12
- Sphinx builder
- MyST Markdown parser
- Math rendering via MathJax (enabled in conf.py)

## Notes

- **Math Support**: Enabled via `sphinx.ext.mathjax` and MyST's `dollarmath` extension
- **Theme**: Using sphinx_rtd_theme for Read the Docs compatibility
- **Format**: Markdown (MyST) for easy editing
- **Images**: Can be added to `_static/` directory and referenced in docs

## Updating Documentation

To add new pages:
1. Create a new `.md` file in `docs/`
2. Add it to the `toctree` in `index.md`
3. Rebuild: `make html`

To modify existing content:
1. Edit the relevant `.md` file
2. Rebuild: `make html`
3. Check `_build/html/` for changes

## Help

For Sphinx documentation:
- https://www.sphinx-doc.org/
- https://myst-parser.readthedocs.io/

For Read the Docs:
- https://docs.readthedocs.io/
