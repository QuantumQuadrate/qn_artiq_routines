"""
I think the feedback does not work well on MOT5 and 6. At some set points (like 0.077V), the RF power does not chnage at all
if I disconnect PD5 from the Sampler, though the measured value on MOT feedback power applet drops to zero!!!
Also, when atom loading is bad, MOT5 and 6 RF powers are clearly off. Also, when I measure MOT5 and 6 powers
at random times by openning the box, they have different values from 200uW to 500uW, totally random. However,
this last measurement might be affected by the RP light not cooling, which is not stabilized.

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

        self.fiber_AOMs = [
            self.dds_AOM_A1,
            self.dds_AOM_A2,
            self.dds_AOM_A3,
            self.dds_AOM_A4,
            self.dds_AOM_A5,
            self.dds_AOM_A6,
        ]
        self.sampler_ch = [7, 5, 3, 4, 6, 0] ### corresponding to how the PDs are connected to the Sampler channel

        self.PD_reading = 0.0

        self.n_average = 5

        self.dt_between_avg_us = 100

        self.Check_power_dataset = ["MOT1_check_power", "MOT2_check_power", "MOT3_check_power",
                                    "MOT4_check_power", "MOT5_check_power", "MOT6_check_power"]

        for i in range(6):
            self.set_dataset(self.Check_power_dataset[i], [0.0], broadcast=True)

        print("prepare - done")

    @kernel
    def _all_fiber_AOMs_off(self):
        for dds in self.fiber_AOMs:
            dds.sw.off()

    @kernel
    def _read_sampler0_avg(self, ch: TInt32, n_avg: TInt32, dt_us: TInt32) -> TFloat:
        SampChVal = 0.0
        buf = np.full(8, 0.0)
        for _ in range(n_avg):
            self.sampler0.sample(buf)
            SampChVal += buf[ch]
            delay(dt_us * us)
        return SampChVal / n_avg

    @kernel
    def run(self):
        self.base.initialize_hardware()
        self.base.initialize_datasets()
        delay(10*ms)

        self.set_dataset("t_step_ms", self.t_step_ms, broadcast=True, persist=False)

        # self.set_point_PD5_AOM_A5 = 0.2

        delay(10 * ms)

        for n in range(self.n_measurements):
            if n % 10 == 0:
                self.print_async(n)
                delay(10*ms)

            ####################################  feedback

            ################## Method 1: using laser_stabilizer.run() to run the feedback on all fiber AOMs.
            self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
            self.laser_stabilizer.run()
            delay(1*ms)


            # ################# Method 2: running stabilizer_AOM_A5, for example, on individual channels
            # ### This works well. Remember, this turns on the object AOM, for example AOM_A5, and leaves it on
            # ### at the end. It does nothing on any other AOM.
            #
            # self.ttl_repump_switch.on()  # Turn off MOT RP
            # self.ttl_pumping_repump_switch.on()  # Turn off Pumping Repump
            #
            # self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
            # delay(1 * ms)
            # self.dds_cooling_DP.sw.on()
            # delay(1 * ms)
            #
            # self.stabilizer_AOM_A1.run(setpoint_index=0)
            # delay(10 * us)
            # self.dds_AOM_A1.sw.off()
            #
            # self.stabilizer_AOM_A2.run(setpoint_index=0)
            # delay(10 * us)
            # self.dds_AOM_A2.sw.off()
            #
            # self.stabilizer_AOM_A3.run(setpoint_index=0)
            # delay(10 * us)
            # self.dds_AOM_A3.sw.off()
            #
            # self.stabilizer_AOM_A4.run(setpoint_index=0)
            # delay(10 * us)
            # self.dds_AOM_A4.sw.off()
            #
            # self.stabilizer_AOM_A5.run(setpoint_index=0)
            # delay(10 * us)
            # self.dds_AOM_A5.sw.off()
            #
            # self.stabilizer_AOM_A6.run(setpoint_index=0)
            # delay(10 * us)
            # self.dds_AOM_A6.sw.off()

            ####################################  Measurement after feedback

            self.dds_cooling_DP.sw.on()
            self.ttl_repump_switch.on()  # Turn off MOT RP
            self.ttl_pumping_repump_switch.on()  # Turn off Pumping Repump
            delay(1 * ms)
            self._all_fiber_AOMs_off()
            delay(1 * ms)

            for i in range(6):
                self.fiber_AOMs[i].sw.on()
                delay(1 * ms)
                self.PD_reading = float(self._read_sampler0_avg(self.sampler_ch[i], int(self.n_average), int(self.dt_between_avg_us)))

                self.append_to_dataset(self.Check_power_dataset[i], self.PD_reading)
                delay(1 * ms)
                self.fiber_AOMs[i].sw.off()
                delay(1 * ms)

            delay(self.t_step_ms * ms)

        ### finally, in case the worker refuses to die
        self.write_results()
        print("*************   Experiment finished   *************")
