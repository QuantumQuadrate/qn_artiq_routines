"""
I think the feedback does not work well on MOT5 and 6. At some set points, the RF power does not chnage at all
if I disconnect PD5 from the Sampler, though the measured value on MOT feedback power applet frops to zero!!!
Also, when atom loading is bad, MOT5 and 6 RF powers are clearly off. Also, when I measure MOT5 and 6 powers
at random times by openning the box, they have different values from 200uW to 500uW, totally random.

So, I wrote this test experiment to run the feedback every couple of seconds, and then measure the PD
values with the sampler and save them in a dataset to plot. The PD values should be stable.

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
from subroutines.experiment_functions import measure_Magnetometer

class Test_feedback_on_MOT5(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.setattr_argument("n_measurements", NumberValue(100, type='int', scale=1, ndecimals=0, step=1))
        self.setattr_argument("t_step_ms", NumberValue(2000, type='int', scale=1, ndecimals=0, step=1))

        self.base.set_datasets_from_gui_args()
        print("build - done")

    def prepare(self):
        """
        performs initial calculations and sets parameter values before
        running the experiment. also sets data filename for now.

        any conversions from human-readable units to machine units (mu) are done here
        """
        self.base.prepare()

        self.n_channels = 8

        self.n_average = 5

        self.smp0 = np.zeros(self.n_channels, dtype=float)
        self.smp1 = np.zeros(self.n_channels, dtype=float)
        self.smp2 = np.zeros(self.n_channels, dtype=float)

        self.avg0 = np.zeros(self.n_channels, dtype=float)
        self.avg1 = np.zeros(self.n_channels, dtype=float)
        self.avg2 = np.zeros(self.n_channels, dtype=float)

        self.smp0List = [self.avg0]
        self.smp1List = [self.avg1]
        self.smp2List = [self.avg2]

        print("prepare - done")

    @kernel
    def run(self):
        self.base.initialize_hardware()
        self.base.initialize_datasets()
        delay(10*ms)

        self.set_dataset("Sampler0Values", self.smp0List, broadcast=True, persist=False)
        self.set_dataset("Sampler1Values", self.smp1List, broadcast=True, persist=False)
        self.set_dataset("Sampler2Values", self.smp2List, broadcast=True, persist=False)
        self.set_dataset("t_step_ms", self.t_step_ms, broadcast=True, persist=False)

        # self.set_point_PD5_AOM_A5 = 0.2

        delay(10 * ms)

        for n in range(self.n_measurements):
            if n % 100 == 0:
                self.print_async(n)
                delay(10*ms)

            self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
            self.laser_stabilizer.run()
            delay(1*ms)
            # self.stabilizer_FORT.run(setpoint_index=1)

            self.dds_cooling_DP.sw.on()  ### Turn on cooling
            delay(1*ms)
            self.dds_AOM_A5.sw.on()
            delay(1*ms)
            self.stabilizer_AOM_A5.run(setpoint_index=0)
            delay(1*ms)

            self.dds_cooling_DP.sw.on()  ### Turn on cooling
            delay(1 * ms)
            self.dds_AOM_A6.sw.on()
            delay(1 * ms)
            self.stabilizer_AOM_A6.run(setpoint_index=0)
            delay(1*ms)

            # self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
            # self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_PGC)

            delay(1*ms)
            self.dds_cooling_DP.sw.on()  ### Turn on cooling

            self.dds_AOM_A5.sw.on()
            self.dds_AOM_A6.sw.on()

            delay(10 * us)

            dummy0 = np.full(8, 0.0)
            # dummy1 = np.full(8, 0.0)
            # dummy2 = np.full(8, 0.0)

            for j in range(self.n_average):
                self.sampler0.sample(self.smp0)  # runs sampler and saves to list
                delay(0.1 * ms)
                # self.sampler1.sample(self.smp1)  # runs sampler and saves to list
                # delay(0.1 * ms)
                # self.sampler2.sample(self.smp2)  # runs sampler and saves to list
                # delay(0.1 * ms)
                dummy0 += self.smp0
                # dummy1 += self.smp1
                # dummy2 += self.smp2
                delay(10 * ms)

            dummy0 /= self.n_average
            # dummy1 /= self.n_average
            # dummy2 /= self.n_average

            self.avg0 = dummy0
            # self.avg1 = dummy1
            # self.avg2 = dummy2

            self.append_to_dataset("Sampler0Values", self.avg0)
            # self.append_to_dataset("Sampler1Values", self.avg1)
            # self.append_to_dataset("Sampler2Values", self.avg2)

            delay(self.t_step_ms * ms)


        print("*************   Experiment finished   *************")
