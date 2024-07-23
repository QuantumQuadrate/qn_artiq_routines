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
from subroutines.experiment_functions import FORT_monitoring_with_Luca_experiment


class MonitorFORTWithLuca(EnvExperiment):

    def build(self):
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.setattr_argument("n_measurements",
                              NumberValue(10, type='int', ndecimals=0, scale=1, step=1))

        self.setattr_argument("t_Luca_exposure", NumberValue(1 * us, ndecimals=3, unit='us'))
        self.setattr_argument("MOT_beams_off", NumberValue(1 * us, ndecimals=3, unit='us'))

        # the Luca will be setup here using pylablib and the analysis will be carried
        # out in self.analyze
        # self.setattr_argument("enable_feedback", BooleanValue(True))

        self.base.set_datasets_from_gui_args()

    def prepare(self):
        self.base.prepare()

    @kernel
    def run(self):
        self.base.initialize_hardware()

        # if not self.control_Luca_in_software:
        print("*****************************  REMEMBER TO START THE CAMERA ACQUISITION  *****************************")
        delay(1000*ms)

        FORT_monitoring_with_Luca_experiment(self)

