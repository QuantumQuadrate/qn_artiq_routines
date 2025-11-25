"""
Testing the chopped RO using self.ttl0._set_sensitivity(0) and (1) in DMA.
Here, I am chopping the FORT and change the ttl sensitivity both on DMA. I am simply measuring the Raman noise from the FORT. The code below
works well and I can do RO when the FORT is off or when the FORT is ON. I have tried cycles up to 4000 which takes about 7ms in total.

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

    @kernel
    def record_chopped_FORT_and_TTL0(self):
        t_FORT_offset = self.core.seconds_to_mu(100 * ns)  ### delay wrt "t0" in each cycle
        t_FORT_OFF = self.core.seconds_to_mu(self.chop_RO_FORT_off) ### the length of the FORT OFF time

        t_cooling_offset = self.core.seconds_to_mu(self.chop_RO_pulse_offset) ### delay wrt "t0" in each cycle
        t_cooling_length = self.core.seconds_to_mu(self.chop_RO_pulse_length) ### duration of the cooling pulse

        t_gate_offset = self.core.seconds_to_mu(self.chop_RO_gate_offset) ### delay wrt t0 for gate to open
        t_gate_length = self.core.seconds_to_mu(self.chop_RO_gate_length) ### gate open time

        with self.core_dma.record("FORT_and_TTL0_chopped"):
            for i in range(4000): ### number of cycles to repeat
                t0 = now_mu()

                at_mu(t0 + t_FORT_offset)
                self.dds_FORT.sw.off()
                at_mu(t0 + t_FORT_offset + t_FORT_OFF)
                self.dds_FORT.sw.on()

                at_mu(t0 + t_cooling_offset)
                self.dds_cooling_DP.sw.on()  ### Turn on cooling
                at_mu(t0 + t_cooling_offset + t_cooling_length)
                self.dds_cooling_DP.sw.off()  ### Turn off cooling

                at_mu(t0 + t_gate_offset)
                # self.ttl_GRIN1_switch.on() ### GRIN1_switch is simply used with set_sensitivity to represent that on scope.
                self.ttl0._set_sensitivity(1)
                self.ttl1._set_sensitivity(1)

                at_mu(t0 + t_gate_offset + t_gate_length)
                # self.ttl_GRIN1_switch.off()
                self.ttl0._set_sensitivity(0)
                self.ttl1._set_sensitivity(0)

                delay(1 * us)

    @kernel
    def run(self):
        self.core.reset()
        self.base.initialize_hardware()

        self.dds_cooling_DP.sw.off()  ### Turn off cooling

        delay(1 * ms)
        self.dds_FORT.sw.on()
        self.dds_AOM_A1.sw.on()
        self.dds_AOM_A2.sw.on()
        self.dds_AOM_A3.sw.on()
        self.dds_AOM_A4.sw.on()
        delay(1 * ms)
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.on()
            self.dds_AOM_A6.sw.on()

        delay(1 * ms)
        self.ttl0._set_sensitivity(0)
        self.ttl_GRIN1_switch.off()

        self.ttl_SPCM0_logic.on() ### turn on the ttl to let the signal pass the logic gates
        self.ttl_SPCM1_logic.on()
        delay(1 * ms)

        ### record and get handle
        self.record_chopped_FORT_and_TTL0()
        delay(10 * ms)
        chopped_handle = self.core_dma.get_handle("FORT_and_TTL0_chopped")
        delay(10 * ms)

        self.ttl_GRIN1_switch.pulse(100*ns)
        delay(5*us)

        ### DMA playback
        self.core_dma.playback_handle(chopped_handle)

        delay(10 * ms)
        SPCM0_counts = self.ttl0.count(now_mu())
        SPCM1_counts = self.ttl1.count(now_mu())

        delay(10 * ms)
        self.dds_FORT.sw.off()
        self.print_async(SPCM0_counts)
        self.print_async(SPCM1_counts)
        print("code done!")
