"""
This code pulses the FORT in the same for our typical single atom experiment cycle,
triggering the Luca camera to take images of the chip to monitor the FORT scattering.

This is for diagnosing whether it will be viable to feedback to the amount of Raman
scattering at 780 nm from the FORT light in the SM fiber that makes it back to the SPCM,
as we want to confirm that that is correlated with the scattering we see in the chamber.

In real experiments we feedback to the FORT in two stages, the latter of which occurs while the MOT is on:
1) the FORT is brought to the level for loading at the beginning of an experiment, when the MOT light is
off and 2) the FORT light is brought down to the science level just before the first atom readout.
"""

from artiq.experiment import *
import logging
import numpy as np
import csv
from time import sleep

import sys, os

cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd + "\\repository\\qn_artiq_routines")

from subroutines.experiment_functions import *


class MonitorFORTWithLuca(EnvExperiment):

    def build(self):
        self.setattr_device('scheduler')
        self.setattr_argument("n_measurements",
                              NumberValue(10, type='int', ndecimals=0, scale=1, step=1))
        self.setattr_argument("warm_up_shots",
                              NumberValue(0, type='int', ndecimals=0, scale=1, step=1))
        self.setattr_argument("MOT_repump_off", BooleanValue(True))
        self.setattr_argument("MOT_light_off", BooleanValue(False))

    def prepare(self):
        # inject variable attributes defined here into the scheduled experiment
        variable_dict = {'MOT_repump_off': self.MOT_repump_off,
                         'MOT_light_off': self.MOT_light_off,
                         'warm_up_shots': self.warm_up_shots,
                         'counts_FORT_loading': 0,
                         'counts_FORT_and_MOT': 0,
                         'counts_FORT_science': 0,
                         'APD_FORT_volts_loading': 0.0,
                         'APD_FORT_volts_science': 0.0,
                         'APD_buffer': np.zeros(8)
                         }

        self.new_job_expid = {'log_level': 30,
                              'file': 'qn_artiq_routines\\GeneralVariableScan.py',
                              'class_name': 'GeneralVariableScan',
                              'arguments': {'n_measurements': str(self.n_measurements),
                                            'scan_variable1_name': 'dummy_variable',
                                            'scan_sequence1': '[1]',
                                            'scan_variable2_name': '',
                                            'scan_sequence2': '',
                                            'override_ExperimentVariables': str(variable_dict),
                                            'experiment_function': 'FORT_monitoring_with_Luca_experiment'},
                              'repo_rev': 'N/A'}

    def run(self):
        self.scheduler.submit(pipeline_name=self.scheduler.pipeline_name,
                              expid=self.new_job_expid,
                              priority=self.scheduler.priority,
                              due_date=None,
                              flush=False)
