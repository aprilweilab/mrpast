import subprocess, os
import os
import sys
sys.path.insert(0, os.path.abspath('../'))

extensions = ["sphinx.ext.autodoc", "sphinx.ext.autosummary"]
html_theme = "pydata_sphinx_theme"
autosummary_generate = True
project = "mrpast"
