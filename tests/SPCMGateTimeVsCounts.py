"""

SPCM Counts vs. Gate Window time


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


class SPCMGateTimeVsCounts(EnvExperiment):

    def build(self):
        self.base = BaseExperiment(experiment=self)
        self.base.build()
        # sets ttls

        self.setattr_argument("n_measurements", NumberValue(1000, ndecimals=0, step=1))
        # self.setattr_argument("t_pulse_mu_list", StringValue('[10,100,1000,10000,100000, 1000000, 10000000]'))
        self.setattr_argument("t_pulse_mu_list", StringValue('[100,200,300,400,500, 600]'))

        self.base.set_datasets_from_gui_args()
        print("build - done")

    def prepare(self):
        self.base.prepare()

        self.t_pulse_mu_list = eval(self.t_pulse_mu_list)
        print(self.t_pulse_mu_list)

    @kernel
    def initialize_hardware(self):      # hardware initialization and setting of ttl switches, and set datasets
        self.base.initialize_hardware()

    def initialize_datasets(self):
        self.base.initialize_datasets()

        # self.set_dataset(self.SPCM0_rate_dataset, [0.0], broadcast=True)

        self.set_dataset("iteration", 0, broadcast=True)
        self.set_dataset("SPCM0_counts_array", [0.0], broadcast=True)
        self.set_dataset("SPCM1_counts_array", [0.0], broadcast=True)


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

        # delay(10*ms)

        # self.dds_cooling_DP.sw.off()
        # self.ttl_repump_switch.on()
        # self.dds_pumping_repump.sw.off()

        delay(10 * ms)

        self.dds_FORT.sw.off()

        self.ttl_excitation_switch.on()
        self.dds_excitation.sw.off()
        self.dds_D1_pumping_DP.sw.off()

        delay(10*ms)


    @kernel
    def experiment_fun(self):
        delay(10 * ms)
        # print(self.t_pulse_mu_list)

        n_cycles = len(self.t_pulse_mu_list)

        SPCM0_counts_array = [0.0] * n_cycles
        SPCM1_counts_array = [0.0] * n_cycles

        self.ttl_SPCM_gate.on() ### remove: SPCM gate no longer exists

        self.ttl_SPCM0._set_sensitivity(1)
        self.ttl_SPCM1._set_sensitivity(1)

        n_step = 0

        for t_pulse_mu in self.t_pulse_mu_list:

            delay(10 * ms)
            print(t_pulse_mu)
            delay(10 * ms)

            SPCM0_counts_sum = 0.0
            SPCM1_counts_sum = 0.0

            for i in range(self.n_measurements):

                now = now_mu()

                at_mu(now)
                self.ttl_SPCM_gate.off()  # turn on ### remove: SPCM gate no longer exists

                at_mu(now + t_pulse_mu)
                self.ttl_SPCM_gate.on()  # turn off ### remove: SPCM gate no longer exists

                delay(1 * us)

                count_now = now_mu()

                SPCM0_counts = self.ttl_SPCM0.count(count_now)
                SPCM1_counts = self.ttl_SPCM1.count(count_now)  # the number of clicks

                SPCM0_counts_sum += SPCM0_counts
                SPCM1_counts_sum += SPCM1_counts

                delay(1*ms)

            SPCM0_counts_array[n_step] = SPCM0_counts_sum / self.n_measurements
            SPCM1_counts_array[n_step] = SPCM1_counts_sum / self.n_measurements

            delay(10 * ms)

            n_step += 1

        self.ttl_SPCM0._set_sensitivity(0)
        self.ttl_SPCM1._set_sensitivity(0)

        for val in SPCM0_counts_array:
            self.append_to_dataset('SPCM0_counts_array', val)
            print(val)

        for val in SPCM1_counts_array:
            self.append_to_dataset('SPCM1_counts_array', val)

        delay(10*ms)


        #     # the ttl.gate_rising(duration) function is equivalent to:
        #     #     ttl._set_sensitivity(1)
        #     #     delay(duration)
        #     #     ttl._set_sensitivity(0)
        #     #     return now_mu()

        # gate_rising
        # self.ttl_SPCM_gate.off() ### remove: SPCM gate no longer exists
        # t_gate_end = self.ttl_SPCM0.gate_rising(self.t_SPCM_second_shot)
        # self.SPCM0_RO2 = self.ttl_SPCM0.count(t_gate_end)


    def run(self):

        self.initialize_hardware()
        self.initialize_datasets()
        self.turn_off_everything()
        self.experiment_fun()





