"""
This code allows for tuning the coils with the homemade potentiometer box
by reading its output into the Sampler and outputting a voltage
from a corresponding Zotino channel. In this way the MOT can be optimized
manually and then the ARTIQ variables for the coil voltages can be updated
automatically.
"""

from artiq.experiment import *
import numpy as np

import sys, os
sys.path.append('C:\\Networking Experiment\\artiq codes\\artiq-master\\repository\\qn_artiq_routines\\')

from utilities.BaseExperiment import BaseExperiment

class MonitorSPCMinApplet(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.setattr_argument("run_time_minutes", NumberValue(1))
        group = "SPCM settings"
         # exposure time of the SPCM
        self.setattr_argument("dt_exposure", NumberValue(300 * ms, unit='ms'), group)
        # saturation limit of the SPCM in counts/s. Can be increased to 10**7 safely, but not higher than 3*10**7.
        self.setattr_argument("sat1s", NumberValue(1 * 10 ** 5), group)  # saturation limit in counts/dt.
        self.setattr_argument("print_count_rate", BooleanValue(True), group)

        self.base.set_datasets_from_gui_args()

        print("build - done")

    def prepare(self):
        """
        performs initial calculations and sets parameter values before
        running the experiment. also sets data filename for now.

        any conversions from human-readable units to machine units (mu) are done here
        """

        self.base.prepare()

        print(self.dt_exposure)
        self.n_steps = int(60*self.run_time_minutes/self.dt_exposure+0.5)
        print(self.n_steps)

        self.count_rate_dataset = 'SPCM_counts_per_s'
        self.set_dataset(self.count_rate_dataset,
                             [0.0],
                             broadcast=True)

        print("prepare - done")

    @kernel
    def run(self):
        self.base.initialize_hardware()

        delay(10 * ms)

        for i in range(self.n_steps):

            tend1 = self.ttl0.gate_rising(self.dt_exposure)
            count1 = self.ttl0.count(tend1)
            count_rate_Hz = count1 / self.dt_exposure
            if self.print_count_rate:
                print(round(count_rate_Hz))
            delay(10 * ms)
            self.append_to_dataset(self.count_rate_dataset, count_rate_Hz)


        print("Experiment finished.")
