# Sphinx configuration for hypar documentation

import os
import sys
from datetime import datetime

# If the package is importable from the repo root, add it to sys.path (adjust if needed)
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, ROOT)

project = "hypar"
author = "hypar contributors"
copyright = f"{datetime.utcnow().year}, {author}"

# Try to import package version if available
try:
    from hypar import __version__ as release  # type: ignore
except Exception:
    release = "0.0.0"

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_autodoc_typehints",
]

myst_enable_extensions = [
    "deflist",
    "colon_fence",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_title = f"{project} docs"

# Autodoc settings
autodoc_member_order = "bysource"
autodoc_typehints = "description"
