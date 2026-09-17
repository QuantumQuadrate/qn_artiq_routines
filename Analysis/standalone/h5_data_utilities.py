"""Forwarding stub so notebooks one directory deeper keep working unchanged.

Every analysis notebook starts with sys.path.append("..") and then
`from h5_data_utilities import *`. Before the standalone/ reorganisation that
".." resolved to Analysis/, where the real module lives. Now that the folders
sit at Analysis/standalone/<folder>/ (and Analysis/standalone/node2/), ".."
resolves to Analysis/standalone/ instead -- this file -- so the import keeps
working and not one notebook had to be edited.

Why it loads the real module by PATH instead of just importing it:

    import sys, os
    sys.path.append(os.path.dirname(os.path.dirname(...)))
    from h5_data_utilities import *          # WRONG

would be a circular self-import. This stub is itself named
h5_data_utilities, so by the time its body runs, Python has already put it in
sys.modules under that name; the inner import would resolve to this
half-initialised module and export nothing. Loading Analysis/
h5_data_utilities.py explicitly under a private module name avoids the
collision entirely.

Keep this file in sync with nothing -- it has no content of its own. The one
real module is Analysis/h5_data_utilities.py.
"""

import importlib.util as _importlib_util
import os as _os

_REAL_PATH = _os.path.join(
    _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))),
    "h5_data_utilities.py",
)

if not _os.path.exists(_REAL_PATH):
    raise ImportError(
        f"h5_data_utilities.py not found at {_REAL_PATH}. This stub expects "
        "the real module one directory up, in Analysis/."
    )

_spec = _importlib_util.spec_from_file_location(
    "_h5_data_utilities_real", _REAL_PATH
)
_real = _importlib_util.module_from_spec(_spec)
_spec.loader.exec_module(_real)

# Re-export everything public, exactly as `import *` from the real module
# would have done.
globals().update(
    {_name: _value for _name, _value in vars(_real).items()
     if not _name.startswith("_")}
)
