"""
This code allows for tuning the balance between each pair of beams with the
homemade potentiometer box by reading its output into the Sampler and
adjusting the set points of a given MOT beam pair. Atom loading can be
optimized in steady state using this routine, in the same way that
SamplerMOTCoilTune can be used.
"""

from artiq.experiment import *
import numpy as np
import logging

import sys, os
# get the current working directory
current_working_directory = os.getcwd()
cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")

from utilities.BaseExperiment import BaseExperiment


class SamplerMOTBeamBalanceTune(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.setattr_argument("FORT_AOM_on", BooleanValue(False))
        self.setattr_argument("set_current_power_setpoints_at_finish", BooleanValue(False))
        self.setattr_argument("leave_coils_on_at_finish", BooleanValue(True))
        self.setattr_argument("run_time_minutes", NumberValue(1))
        self.setattr_argument("max_setpoint_percent_deviation",
                              NumberValue(0.1))

        group = "SPCM settings"
         # exposure time of the SPCM
        self.setattr_argument("dt_exposure", NumberValue(15 * ms, unit='ms'), group)

        # when to run the AOM feedback (after how many iterations in the for loops)
        self.setattr_argument("AOM_feedback_period_cycles", NumberValue(500), "Laser feedback")
        self.setattr_argument("monitor_only", BooleanValue(False), "Laser feedback")

        self.base.set_datasets_from_gui_args()
        print("build - done")

    def prepare(self):
        """
        performs initial calculations and sets parameter values before
        running the experiment. also sets data filename for now.

        any conversions from human-readable units to machine units (mu) are done here
        """

        self.base.prepare()

        self.n_steps = int(60*self.run_time_minutes/self.dt_exposure+0.5)
        self.sampler_buffer = np.zeros(8)
        self.sampler_channels = [0,1,2] # the sampler channels to read

        self.setpoint_datasets = ["set_point_PD1_AOM_A1","set_point_PD2_AOM_A2","set_point_PD3_AOM_A3",
                                  "set_point_PD4_AOM_A4","set_point_PD5_AOM_A5","set_point_PD6_AOM_A6"]
        self.default_setpoints = [getattr(self,dataset) for dataset in self.setpoint_datasets]

        self.set_dataset(self.count_rate_dataset,
                         [0.0],
                         broadcast=True)
        print("prepare - done")

    @kernel
    def run(self):
        self.base.initialize_hardware()

        # Turn on AOMs to load the MOT.
        self.dds_cooling_DP.sw.on()
        self.dds_AOM_A2.sw.on()
        self.dds_AOM_A3.sw.on()
        self.dds_AOM_A1.sw.on()
        self.dds_AOM_A6.sw.on()
        self.dds_AOM_A4.sw.on()
        self.dds_AOM_A5.sw.on()

        if self.FORT_AOM_on:
            self.dds_FORT.sw.on()
        else:
            self.dds_FORT.sw.off()

            # wait for AOMs to thermalize
            delay(3000 * ms)

        delay(1 * ms)

        self.laser_stabilizer.run(monitor_only=self.monitor_only)
        if self.FORT_AOM_on:
            self.dds_FORT.sw.on()

        print("ready!")

        for i in range(self.n_steps):

            if (i % self.AOM_feedback_period_cycles) == 0:
                print("running feedback")
                self.core.break_realtime()
                self.laser_stabilizer.run(monitor_only=self.monitor_only)
                delay(100*ms)
                if self.FORT_AOM_on:
                    self.dds_FORT.sw.on()
            else:
                delay(10*ms)

            t_end = self.ttl0.gate_rising(self.dt_exposure)
            counts = self.ttl0.count(t_end)
            count_rate_per_s = counts / self.dt_exposure

            delay(1 * ms)
            self.append_to_dataset(self.count_rate_dataset, count_rate_per_s)

            self.sampler1.sample(self.sampler_buffer)
            sp1_factor = 1 + self.sampler_buffer[self.sampler_channels[0]] * self.max_setpoint_percent_deviation/3.5
            sp2_factor = 1 - self.sampler_buffer[self.sampler_channels[0]] * self.max_setpoint_percent_deviation/3.5
            sp3_factor = 1 + self.sampler_buffer[self.sampler_channels[1]] * self.max_setpoint_percent_deviation/3.5
            sp4_factor = 1 - self.sampler_buffer[self.sampler_channels[1]] * self.max_setpoint_percent_deviation/3.5
            sp5_factor = 1 + self.sampler_buffer[self.sampler_channels[2]] * self.max_setpoint_percent_deviation/3.5
            sp6_factor = 1 - self.sampler_buffer[self.sampler_channels[2]] * self.max_setpoint_percent_deviation/3.5
            self.print_async(sp1_factor,sp2_factor,sp3_factor,sp4_factor,sp5_factor,sp6_factor)

            # update the fiber AOM setpoints
            self.stabilizer_AOM_A1.set_point = self.default_setpoints[0]*sp1_factor
            self.stabilizer_AOM_A2.set_point = self.default_setpoints[1]*sp2_factor
            self.stabilizer_AOM_A3.set_point = self.default_setpoints[2]*sp3_factor
            self.stabilizer_AOM_A4.set_point = self.default_setpoints[3]*sp4_factor
            self.stabilizer_AOM_A5.set_point = self.default_setpoints[4]*sp5_factor
            self.stabilizer_AOM_A6.set_point = self.default_setpoints[5]*sp6_factor

            delay(1*ms)

        if not self.leave_coils_on_at_finish:
            for ch in self.coil_channels:
                self.zotino0.write_dac(ch, 0.0)
                self.zotino0.load()
                delay(1 * ms)

        if self.set_current_power_setpoints_at_finish:
            current_power_setpoints = [
                self.stabilizer_AOM_A1.set_point,
                self.stabilizer_AOM_A2.set_point,
                self.stabilizer_AOM_A3.set_point,
                self.stabilizer_AOM_A4.set_point,
                self.stabilizer_AOM_A5.set_point,
                self.stabilizer_AOM_A6.set_point
            ]

            for i in range(6):
                self.set_dataset(self.setpoint_datasets[i], current_power_setpoints[i], broadcast=True, persist=True)

        print("Experiment finished.")
