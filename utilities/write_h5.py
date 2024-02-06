"""
for writing an h5 file from within an ARTIQ experiment, i.e., before the worker has exited.
"""

import sys
import time
import os
import inspect
import logging
import traceback
from collections import OrderedDict
import importlib.util
import linecache

import h5py

from sipyco import pipe_ipc, pyon
from sipyco.packed_exceptions import raise_packed_exc
from sipyco.logging_tools import multiline_log_config

from artiq.experiment import *
from artiq import tools
from artiq.master.worker_db import DeviceManager, DatasetManager
from artiq.language.environment import (
    is_public_experiment, TraceArgumentManager, ProcessArgumentManager
)
from artiq.language.core import set_watchdog_factory, TerminationRequested
from artiq.language.types import TBool
from artiq.compiler import import_cache
from artiq.coredevice.core import CompileError, host_only, _render_diagnostic
from artiq import __version__ as artiq_version 


@rpc(flags={"async"})
def write_results(experiment, name=None):
    try:    
        experiment.setattr_device("scheduler")
        rid = experiment.scheduler.rid
        run_time = 0
        ##Tracking time of last update
        start_time = time.time()
        expid = experiment.scheduler.expid
        dataset_mgr = experiment._HasEnvironment__dataset_mgr
        ## rid_counter.py in master relies on the previously used RID if the RID cache is cleared.
        if name != None:
            filename = "{:09}-{}_{}.h5".format(rid, experiment.__class__.__name__, name)
        else:
            filename = "{:09}-{}Synchronous.h5".format(rid, experiment.__class__.__name__)
        with h5py.File(filename, "w") as f:
            dataset_mgr.write_hdf5(f)
            f["artiq_version"] = artiq_version
            f["rid"] = rid
            f["start_time"] = start_time
            f["run_time"] = run_time
            f["expid"] = pyon.encode(expid)
    except Exception as e:
        pass
    