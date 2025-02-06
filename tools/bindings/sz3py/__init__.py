import os
import sys
import importlib.util

# Get the directory of this package
_package_dir = os.path.dirname(__file__)

# Find the .so file
_so_file = [f for f in os.listdir(_package_dir) if f.endswith('.so')]
if not _so_file:
    raise ImportError("Could not find the compiled sz3py shared library (.so)")

# Load the .so file dynamically
_so_path = os.path.join(_package_dir, _so_file[0])
_spec = importlib.util.spec_from_file_location("sz3py", _so_path)
_sz3py = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_sz3py)

# Import everything from the .so module
globals().update({name: getattr(_sz3py, name) for name in dir(_sz3py) if not name.startswith("__")})
