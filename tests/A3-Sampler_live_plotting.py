"""
Monitor three Sampler channels for live plotting in applet, using Sampler1 only for now.
Choose the channel below (hard coded).

"""

from artiq.experiment import *
import numpy as np
from datetime import datetime as dt

import sys, os
### get the current working directory
current_working_directory = os.getcwd()
cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")

from utilities.BaseExperiment import BaseExperiment

class Sampler_live_plotting(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.setattr_argument("n_measurements", NumberValue(100, type='int', scale=1, ndecimals=0, step=1))
        self.setattr_argument("n_average", NumberValue(10, type='int', scale=1, ndecimals=0, step=1))
        self.setattr_argument("t_step_ms", NumberValue(10, type='int', scale=1, ndecimals=0, step=1))

        self.base.set_datasets_from_gui_args()
        print("build - done")

    def prepare(self):
        """
        performs initial calculations and sets parameter values before
        running the experiment. also sets data filename for now.

        any conversions from human-readable units to machine units (mu) are done here
        """
        self.base.prepare()
        self.Sampler_ch_A = 0
        self.Sampler_ch_B = 1
        self.Sampler_ch_C = 2

        print("prepare - done")

    @kernel
    def run(self):
        self.base.initialize_hardware()
        # self.base.initialize_datasets()
        self.expt()
        print("*************   Experiment finished   *************")

    @kernel
    def expt(self):

        # measure_Magnetometer(self) ### works, but I want to have more control

        ### x,y, and z axes are connected to Sampler2 Ch2,3, and 4, respectively.

        self.set_dataset("Sampler_ch_A", [0.0], broadcast=True)
        self.set_dataset("Sampler_ch_B", [0.0], broadcast=True)
        self.set_dataset("Sampler_ch_C", [0.0], broadcast=True)

        self.core.break_realtime()

        for n_measrement in range(self.n_measurements):
            if n_measrement % max(1, self.n_measurements// 10) == 0:
                self.print_async("progress (%): ", (n_measrement / self.n_measurements) * 100)

            measurement_buf = np.array([0.0] * 8)
            Sampler_ch_A_value = 0.0
            Sampler_ch_B_value = 0.0
            Sampler_ch_C_value = 0.0

            for i in range(self.n_average):
                self.sampler1.sample(measurement_buf)
                Sampler_ch_A_value += measurement_buf[self.Sampler_ch_A]
                Sampler_ch_B_value += measurement_buf[self.Sampler_ch_B]
                Sampler_ch_C_value += measurement_buf[self.Sampler_ch_C]
                delay(100 * us)

            Sampler_ch_A_value /= self.n_average
            Sampler_ch_B_value /= self.n_average
            Sampler_ch_C_value /= self.n_average

            self.append_to_dataset("Sampler_ch_A", Sampler_ch_A_value)
            self.append_to_dataset("Sampler_ch_B", Sampler_ch_B_value)
            self.append_to_dataset("Sampler_ch_C", Sampler_ch_C_value)
            delay(self.t_step_ms * ms)

        ### finally, in case the worker refuses to die
        self.write_results()