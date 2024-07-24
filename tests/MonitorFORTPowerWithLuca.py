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
sys.path.append(cwd+"\\repository\\qn_artiq_routines")

from utilities.BaseExperiment import BaseExperiment
from subroutines.experiment_functions import *


class MonitorFORTWithLuca(EnvExperiment):

    def build(self):
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.setattr_argument("n_measurements",
                              NumberValue(10, type='int', ndecimals=0, scale=1, step=1))

        # configure the Luca in Andor Solis
        self.setattr_argument("t_Luca_exposure", NumberValue(1 * us, ndecimals=3, unit='us'))

        self.base.set_datasets_from_gui_args()

    def prepare(self):
        self.new_job_expid = {'log_level': 30,
                              'file': 'qn_artiq_routines\\GeneralVariableScan.py',
                              'class_name': 'GeneralVariableScan',
                              'arguments': {'n_measurements': str(self.n_measurements),
                                            'scan_variable1': 'dummy_variable',
                                            'scan_sequence1': '[1]',
                                            'scan_variable2': '',
                                            'scan_sequence2': '',
                                            'override_ExperimentVariables':
                                                '{"t_Luca_exposure":'+str(self.t_Luca_exposure)+'}',
                                            'experiment_function': 'FORT_monitoring_with_Luca_experiment'},
                              'repo_rev': 'N/A'}

    def initialize_hardware(self):
        self.base.initialize_hardware()

    def run(self):
        self.scheduler.submit(pipeline_name=self.scheduler.pipeline_name,
                              expid=self.new_job_expid,
                              priority=self.scheduler.priority,
                              due_date=None,
                              flush=False)
