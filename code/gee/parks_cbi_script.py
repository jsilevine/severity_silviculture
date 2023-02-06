#!/usr/bin/env python3

from geemap.conversion import *

js_dir = "./"

# Convert all Earth Engine JavaScripts in a folder recursively to Python scripts.
js_to_python_dir(in_dir=js_dir, out_dir=js_dir, use_qgis=True)
