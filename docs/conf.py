import os
import sys
sys.path.insert(0, os.path.abspath('../'))

extensions = ["sphinx.ext.autodoc", "sphinx.ext.autosummary"]
html_theme = "pydata_sphinx_theme"
autosummary_generate = True
project = "mrpast"

html_theme_options = {
    "external_links": [
        {
            "url": "https://aprilweilab.github.io",
            "name": "Wei Lab Website",
        }
    ],
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/aprilweilab/mrpast",
            "icon": "fa-brands fa-github",
        },
        {
            "name": "PyPI",
            "url": "https://pypi.org/project/mrpast",
            "icon": "fa-custom fa-pypi",
        }
    ],
}
html_static_path = ["_static"]
html_js_files = [
    ("custom-icons.js", {"defer": "defer"}),
]
