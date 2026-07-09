"""Sphinx configuration for pypkatool.

Building these docs does not require the pypka/pKAI conda environment
(environment.yml) - only docs/requirements.txt. `pypkatool.core` checks for
the `pkai` package at import time (see `_require_pkai()`), so a minimal stub
package under `_stubs/` is put on `sys.path` below purely to satisfy that
check; the real `pypka`/`pkai` imports inside the module are all function-
local and are never executed just by importing the module for autodoc.
"""
import sys
from pathlib import Path

DOCS_DIR = Path(__file__).parent
sys.path.insert(0, str(DOCS_DIR.parent))          # so `import pypkatool` resolves
sys.path.insert(0, str(DOCS_DIR / "_stubs"))       # stub `pkai` package for _require_pkai()

project = "pypkatool"
copyright = "2026, Ropón-Palacios G."
author = "Ropón-Palacios G."
release = "1.0.0"

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
]

myst_enable_extensions = ["colon_fence", "deflist"]
source_suffix = {".md": "markdown", ".rst": "restructuredtext"}

autodoc_member_order = "bysource"
autodoc_typehints = "description"

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "_stubs"]

html_theme = "furo"
html_title = "pypkatool"

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
}
