"""

SPCM Counts vs. Gate Window time

Problems:

"""

from artiq.experiment import *
import logging
import numpy as np
import pyvisa as visa

import sys, os
# get the current working directory
cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")

from utilities.BaseExperiment import BaseExperiment
import subroutines.experiment_functions as exp_functions
from subroutines.experiment_functions import load_MOT_and_FORT

from utilities.conversions import dB_to_V_kernel as dB_to_V


class ExcitationSPCMGateStartTime(EnvExperiment):

    def build(self):
        self.base = BaseExperiment(experiment=self)
        self.base.build()
        # sets ttls

        self.setattr_argument("n_measurements", NumberValue(10000000, ndecimals=0, step=1))
        # self.setattr_argument("exc_RF_in_dBm", NumberValue(1., ndecimals=0, step=1))
        self.setattr_argument("exc_pulse_length_mu", NumberValue(20, ndecimals=0, step=1))
        self.setattr_argument("gate_pulse_length_mu", NumberValue(100, ndecimals=0, step=1))


        self.setattr_argument("t_pulse_mu_list", StringValue('[100,200,300,400,500, 600]'))

        self.setattr_argument("test_with_specific_ao_setting", BooleanValue(False))

        self.base.set_datasets_from_gui_args()
        print("build - done")

    def prepare(self):
        self.base.prepare()

        self.t_pulse_mu_list = eval(self.t_pulse_mu_list)

    @kernel
    def initialize_hardware(self):      # hardware initialization and setting of ttl switches, and set datasets
        self.base.initialize_hardware()

    def initialize_datasets(self):
        self.base.initialize_datasets()

        # self.set_dataset(self.count_rate_dataset, [0.0], broadcast=True)

        self.set_dataset("iteration", 0, broadcast=True)
        self.set_dataset("SPCM_counts_array", [0.0], broadcast=True)
        self.set_dataset("SPCM_counts1_array", [0.0], broadcast=True)


    @kernel
    def turn_off_everything(self):

        delay(1*s)

        self.dds_cooling_DP.sw.off()
        delay(1 * ms)

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        delay(1*ms)

        self.dds_cooling_DP.sw.off()
        self.ttl_exc0_switch.on()
        self.dds_pumping_repump.sw.off()

        delay(1 * ms)

        self.dds_FORT.sw.off()

        self.ttl_GRIN1_switch.on()
        self.dds_excitation.sw.off()
        self.dds_D1_pumping_DP.sw.off()

        delay(10*ms)

    @kernel
    def specific_ao_setting(self):

        delay(1*s)
        # Repump
        self.ttl_repump_switch.off()

        # Cooling
        self.dds_cooling_DP.sw.off()

        # Pumping Repump
        self.dds_pumping_repump.sw.off()

        # FORT
        self.dds_FORT.sw.off()

        # Excitation
        self.ttl_GRIN1_switch.on()
        self.dds_excitation.sw.off()

        # D1
        self.dds_D1_pumping_DP.sw.off()

        delay(1 * ms)

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        delay(10*ms)


    @kernel
    def experiment_fun(self):

        self.core.reset()

        delay(10 * ms)

        #turning SPCM gate OFF & Sensitivity ON
        self.ttl_SPCM_gate.on()
        self.ttl_SPCM0._set_sensitivity(1)
        self.ttl_SPCM1._set_sensitivity(1)

        # setting excitation aom to 1dBm
        # self.dds_excitation.set(frequency=self.f_excitation,amplitude=dB_to_V(self.exc_RF_in_dBm))
        self.dds_excitation.set(frequency=self.f_excitation, amplitude=dB_to_V(1.0))
        # don't have to set it back; any code that runs feedback will do

        #turning on excitation; ttl switch still blocked
        self.ttl_exc0_switch.off()
        self.dds_excitation.sw.on()

        delay(10*ms)

        for t_pulse_mu in self.t_pulse_mu_list:

            delay(10 * ms)
            # print(t_pulse_mu)
            delay(10 * ms)

            SPCM_counts_sum = 0.0
            SPCM_counts1_sum = 0.0

            t_exc_on_mu = 500 # 500ns delay before turning excitation on.

            for i in range(self.n_measurements):

                now = now_mu()

                # excitation on
                at_mu(now + t_exc_on_mu)
                self.ttl_GRIN1_switch.off()

                # gate on
                at_mu(now + t_exc_on_mu + t_pulse_mu)
                self.ttl_SPCM_gate.off()  # turn on

                # excitation off
                at_mu(now + t_exc_on_mu + self.exc_pulse_length_mu)
                self.ttl_GRIN1_switch.on()

                # gate off
                at_mu(now + t_exc_on_mu + t_pulse_mu + self.gate_pulse_length_mu)
                self.ttl_SPCM_gate.on()  # turn off

                count_now = now_mu()

                SPCM_counts = self.ttl_SPCM0.count(count_now)
                SPCM_counts1 = self.ttl_SPCM1.count(count_now)  # the number of clicks

                SPCM_counts_sum += SPCM_counts
                SPCM_counts1_sum += SPCM_counts1

                delay(10 * us)

            self.append_to_dataset('SPCM_counts_array', SPCM_counts_sum)
            self.append_to_dataset('SPCM_counts1_array', SPCM_counts1_sum)

            delay(1 * ms)


        self.ttl_SPCM0._set_sensitivity(0)
        self.ttl_SPCM1._set_sensitivity(0)

        # for val in SPCM_counts_array:
        #     self.append_to_dataset('SPCM_counts_array', val)
        #
        # for val in SPCM_counts1_array:
        #     self.append_to_dataset('SPCM_counts1_array', val)

        delay(10*ms)


        #     # the ttl.gate_rising(duration) function is equivalent to:
        #     #     ttl._set_sensitivity(1)
        #     #     delay(duration)
        #     #     ttl._set_sensitivity(0)
        #     #     return now_mu()

        # gate_rising
        # self.ttl_SPCM_gate.off()
        # t_gate_end = self.ttl_SPCM0.gate_rising(self.t_SPCM_second_shot)
        # self.counts2 = self.ttl_SPCM0.count(t_gate_end)


    def run(self):

        self.initialize_hardware()
        self.initialize_datasets()

        if self.test_with_specific_ao_setting:
            self.specific_ao_setting()
        else:
            self.turn_off_everything()

        self.experiment_fun()