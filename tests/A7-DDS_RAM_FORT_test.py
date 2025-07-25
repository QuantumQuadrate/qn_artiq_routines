"""
This is another test experiment that is similar to our actual experiment that uses BaseExperiment and aom_feedback.
The goal is to test that we can run FORT on RAM mode (turning ON/OFF smoothly) without interfering with other AOMs.
This is one level more advance that DDS_RAM_amplitude_test.py.

It works well; we exit the RAM mode to run the stabilizer and re-enter the RAM mode.

Akbar 2025-06-24

"""

from artiq.experiment import *
from artiq.coredevice.ad9910 import RAM_MODE_RAMPUP, RAM_DEST_ASF
from artiq.coredevice.urukul import CFG_MASK_NU
from artiq.language.types import TInt32
from artiq.language import ms, us, ns, MHz
import numpy as np
import math

import sys, os
### get the current working directory
current_working_directory = os.getcwd()
cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")

from utilities.BaseExperiment import BaseExperiment

class DDS_RAM_FORT_test(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.setattr_argument("n_measurements", NumberValue(100, type='int', scale=1, ndecimals=0, step=1))
        self.setattr_argument("n_average", NumberValue(10, type='int', scale=1, ndecimals=0, step=1))

        self.base.set_datasets_from_gui_args()
        print("build - done")

    def prepare(self):
        """
        performs initial calculations and sets parameter values before
        running the experiment. also sets data filename for now.

        any conversions from human-readable units to machine units (mu) are done here
        """
        self.base.prepare()

        print("prepare - done")

    @kernel
    def run(self):
        self.base.initialize_hardware()
        self.base.initialize_datasets()
        self.core.break_realtime()

        if self.enable_laser_feedback:
            ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
            self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
            self.stabilizer_FORT.run(setpoint_index=1)  # the science setpoint
            self.laser_stabilizer.run()


        ################### preparing and writing the RAM
        self.prepare_FORT_RAM_profile()

        self.dds_FORT.set_cfr1(ram_enable=0)  ### disable RAM mode to write the config
        self.dds_FORT.cpld.io_update.pulse_mu(8)  ### pulse the ttl to update and implement settings

        ### Configures the RAM playback engine
        self.dds_FORT.set_profile_ram(
            start=0,
            end=self.FORT_total_points - 1,
            step=self.FORT_step_size,
            profile=7,
            mode=RAM_MODE_RAMPUP,
        )

        ### write the data onto RAM
        self.dds_FORT.cpld.io_update.pulse_mu(8)
        self.dds_FORT.write_ram(self.FORT_amplitudes_list)
        ################### End of preparing and writing the RAM


        ###########  Test 1: turning FORT ON/OFF outside RAM.
        #### note that we HAVE TO set the freq. and ampl. before turning on the AOM.
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
        self.dds_cooling_DP.sw.off()
        self.dds_FORT.sw.on()
        delay(1*ms)
        self.dds_cooling_DP.sw.on()
        self.dds_FORT.sw.off()
        delay(1*ms)
        ###########  End of Test 1


        self.expt()
        print("*************   Experiment finished   *************")

    @kernel
    def expt(self):
        delay(1 * ms)

        ############################################  Turning on FORT in RAM mode
        ### Enabling RAM playback, not playing yet.
        self.dds_FORT.set_cfr1(ram_enable=1,
                               ram_destination=RAM_DEST_ASF)

        ### turns on the FORT to loading set point
        self.dds_FORT.set_profile_ram(
            start=self.FORT_total_points // 2,
            end=self.FORT_total_points - 1,
            step=self.FORT_step_size,
            profile=7,
            mode=RAM_MODE_RAMPUP)
        self.dds_FORT.cpld.io_update.pulse_mu(8)
        self.dds_FORT.sw.on()
        self.dds_FORT.cpld.io_update.pulse_mu(8)
        #######################################  End of Turning on FORT
        delay(2 * ms)

        for measurement in range(self.n_measurements):
            self.ttl7.on()

            ### Configure the RAM to playback the first half (ramping down to science point)
            self.dds_FORT.set_profile_ram(
                start=0,
                end=len(self.FORT_amplitudes_list) // 2 - 1,
                step=self.FORT_step_size,
                profile=7,
                mode=RAM_MODE_RAMPUP)
            self.dds_FORT.cpld.io_update.pulse_mu(8)
            delay(5*ms)

            self.ttl7.off()


            ########## test 2
            ### changing cooling DP AOM set point without affecting FORT:
            self.dds_FORT.cpld.cfg_write(self.dds_FORT.cpld.cfg_reg | 1 << CFG_MASK_NU + 0)  # Mask the DDS channel which is on RAM mode
            self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO,
                                    amplitude=self.ampl_cooling_DP_MOT * self.p_cooling_DP_RO)
            self.dds_FORT.cpld.cfg_write(self.dds_FORT.cpld.cfg_reg & ~(1 << CFG_MASK_NU + 0))  # Unmask the DDS channel which is on RAM mode
            self.dds_cooling_DP.sw.off()
            delay(1*ms)
            self.dds_cooling_DP.sw.on()
            delay(1 * ms)
            ########## End of test 2


            # ########## test 3: exiting and re-entering the RAM mode
            # ### exiting the RAM mode
            # self.dds_FORT.set_cfr1(ram_enable=0)
            # self.dds_FORT.cpld.io_update.pulse_mu(8)
            # delay(1*ms)
            #
            # if self.enable_laser_feedback:
            #     ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
            #     self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
            #     self.stabilizer_FORT.run(setpoint_index=1)  # the science setpoint
            #     self.laser_stabilizer.run()
            #
            # ### Enabling RAM playback again, not playing yet.
            # self.dds_FORT.set_cfr1(ram_enable=1,
            #                        ram_destination=RAM_DEST_ASF)
            # ########## End of test 3
            # ### Disable test 3 to avoid FORT turning off on the scope.


            ### turns on the FORT to loading set point
            self.dds_FORT.set_profile_ram(
                start=self.FORT_total_points // 2,
                end=self.FORT_total_points - 1,
                step=self.FORT_step_size,
                profile=7,
                mode=RAM_MODE_RAMPUP)
            self.dds_FORT.cpld.io_update.pulse_mu(8)
            self.dds_FORT.sw.on()
            self.dds_FORT.cpld.io_update.pulse_mu(8)
            #######################################  End of Turning on FORT

            delay(5*ms)



        ### exiting the RAM mode
        self.dds_FORT.set_cfr1(ram_enable=0)
        self.dds_FORT.cpld.io_update.pulse_mu(8)

    @kernel
    def prepare_FORT_RAM_profile(self):
        """
        Prepares FORT RAM waveform to be used to ramp down and ramp up the FORT smoothly.
        """

        FORT_steps_rise = 50  # number of amplitude points during the rise and fall time
        FORT_amp_high = self.stabilizer_FORT.amplitudes[0]  # low and high amplitudes in scale from 0 to 1
        FORT_amp_low = self.stabilizer_FORT.amplitudes[1]

        # print(FORT_amp_high)
        # delay(10*ms)
        # print(FORT_amp_low)
        # delay(10 * ms)

        FORT_step_ticks = int(
            (self.t_FORT_ramp / FORT_steps_rise) / (
                    4 * ns))  ### Step size

        FORT_dwell_time = 100 * us

        FORT_steps_dwell = int(
            FORT_steps_rise / self.t_FORT_ramp * FORT_dwell_time)  # Number of dwell points

        FORT_x_vals = [i / (FORT_steps_rise - 1) for i in range(FORT_steps_rise)]
        FORT_norm = [6 * x ** 5 - 15 * x ** 4 + 10 * x ** 3 for x in FORT_x_vals]  # quintic smoothstep function

        FORT_amp_points_rise = [FORT_amp_low + n * (FORT_amp_high - FORT_amp_low) for n in FORT_norm]
        FORT_amp_points_fall = [FORT_amp_points_rise[i] for i in range(FORT_steps_rise - 1, -1, -1)]

        ### The full waveform.
        FORT_amp_points = (
                FORT_amp_points_fall +
                [FORT_amp_low] * FORT_steps_dwell +
                FORT_amp_points_rise
        )

        ### some data conversion needed for RAM
        FORT_amplitudes_arr = [0] * len(FORT_amp_points)
        self.dds_FORT.amplitude_to_ram(FORT_amp_points, FORT_amplitudes_arr)
        FORT_amplitudes_list = list(FORT_amplitudes_arr)

        ### This is calculation of steps based on above parameters
        FORT_total_points = len(FORT_amplitudes_list)

        self.FORT_step_size = FORT_step_ticks
        self.FORT_total_points = FORT_total_points
        self.FORT_amplitudes_list = FORT_amplitudes_list

        delay(1 * ms)
