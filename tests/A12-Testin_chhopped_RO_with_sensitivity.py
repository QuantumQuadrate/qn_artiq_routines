"""
Testing the chopped RO using self.ttl0._set_sensitivity(0) and (1) in DMA.
Here, I am chopping the FORT and change the ttl sensitivity both on DMA.

Akbar 2025-11-24

"""

from artiq.experiment import *
import logging
import numpy as np

import sys, os
# get the current working directory
cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")

from utilities.BaseExperiment import BaseExperiment
from utilities.conversions import dB_to_V_kernel as dB_to_V


class Testing_chopped_RO_with_sensitivity(EnvExperiment):

    def build(self):
        self.base = BaseExperiment(experiment=self)
        self.base.build()

    def prepare(self):
        self.base.prepare()

    # @kernel
    # def record_pulses_FORT(self):
    #     with self.core_dma.record("FORT_pulses"):
    #         for i in range(400):
    #             self.dds_FORT.sw.on()
    #             delay(1 * us)
    #             self.dds_FORT.sw.off()
    #             delay(1 * us)
    #
    # @kernel
    # def record_TTL0_pulse_sensitivity(self):
    #     with self.core_dma.record("TTL0_sensivity_pulses"):
    #         for i in range(400):
    #             self.ttl0._set_sensitivity(1)
    #             delay(1 * us)
    #             self.ttl0._set_sensitivity(0)
    #             delay(1 * us)

    @kernel
    def record_chopped_FORT_and_TTL0(self):
        with self.core_dma.record("FORT_and_TTL0_chopped"):
            for i in range(400):
                ### gate on + FORT on
                self.ttl0._set_sensitivity(1)
                self.dds_FORT.sw.on()
                delay(1 * us)

                # gate off + FORT off
                self.ttl0._set_sensitivity(0)
                self.dds_FORT.sw.off()
                delay(1 * us)

    @kernel
    def run(self):
        self.core.reset()
        self.base.initialize_hardware()

        delay(1 * ms)
        self.dds_FORT.sw.on()
        delay(1 * ms)
        self.ttl0._set_sensitivity(0)

        self.ttl_SPCM0_logic.on() ### turn on the ttl to let the signal pass the logic gates
        self.ttl_SPCM1_logic.on()
        delay(1 * ms)

        ### record and get handle
        self.record_chopped_FORT_and_TTL0()
        delay(10 * ms)
        chopped_handle = self.core_dma.get_handle("FORT_and_TTL0_chopped")
        delay(10 * ms)

        ### just one DMA playback, no parallel
        self.core_dma.playback_handle(chopped_handle)

        delay(10 * ms)
        SPCM0_SinglePhoton = self.ttl0.count(now_mu())

        delay(10 * ms)
        self.dds_FORT.sw.off()
        self.print_async(SPCM0_SinglePhoton)


    # @kernel
    # def run(self):
    #     self.core.reset()
    #
    #     self.base.initialize_hardware()
    #
    #     delay(1*ms)
    #     self.dds_FORT.sw.on()  ### turns FORT on
    #     delay(1 * ms)
    #     self.ttl0._set_sensitivity(0)
    #
    #     self.ttl_SPCM0_logic.on() ### turn on the ttl to let the signal pass the logic gates
    #     self.ttl_SPCM1_logic.on()
    #
    #     delay(1 * ms)
    #
    #
    #     #### The following seems to work as we want; sets the sensitivity using DMA and counts the events afterwards.
    #     self.record_TTL0_pulse_sensitivity()
    #     delay(10 * ms)
    #     self.record_pulses_FORT()
    #     delay(10*ms)
    #     TTL0_gates_handle = self.core_dma.get_handle("TTL0_sensivity_pulses")
    #     delay(10 * ms)
    #     FORT_pulses_handle = self.core_dma.get_handle("FORT_pulses")
    #     delay(10 * ms)
    #
    #     with parallel:
    #         self.core_dma.playback_handle(FORT_pulses_handle)
    #         # self.core_dma.playback_handle(TTL0_gates_handle)
    #     delay(10 * ms)
    #     SPCM0_SinglePhoton = self.ttl0.count(now_mu())
    #
    #     delay(10 * ms)
    #     self.dds_FORT.sw.off()  ### turns FORT off
    #     self.print_async(SPCM0_SinglePhoton)
    #
    #
    #     print("code done!")
