"""
Connect ttl_SPCM_gate to an oscilloscope with a tee, sending the other port to an SPCM.
Monitor the SPCM output on the scope too.
"""

from artiq.experiment import *

import sys, os
# get the current working directory
cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")

from utilities.BaseExperiment import BaseExperiment
from subroutines.experiment_functions import load_MOT_and_FORT


class SPCMGateRiseTimeTest(EnvExperiment):

    def build(self):
        self.base = BaseExperiment(experiment=self)
        self.base.build()
        self.setattr_argument("pulse_duty_cycle", NumberValue(0.2, precision=1, step=1))
        self.setattr_argument("pulse_period", NumberValue(500*ns, unit='ns', precision=1, step=1))
        self.setattr_argument("disable_pulse_sequence", BooleanValue(False)) # for viewing SPCM output in steady state

        self.base.set_datasets_from_gui_args()

    def prepare(self):
        self.base.prepare()

    @kernel
    def record_dma_pulses(self):
        period_mu = int(self.pulse_period / ns)
        t_pulse = self.pulse_duty_cycle * self.pulse_period

        with self.core_dma.record("fast_SPCM_gating"):

            start_mu = now_mu()
            for i in range(2000):
                at_mu(start_mu + i * period_mu)
                self.ttl_SPCM_gate.pulse(t_pulse)

    @kernel
    def run(self):
        self.base.initialize_hardware()
        self.record_dma_pulses()

        dma_handle = self.core_dma.get_handle("fast_SPCM_gating")


        self.laser_stabilizer.run()
        delay(1 * ms)

        load_MOT_and_FORT(self)  # we'll measure scattered photons from all the beams, Raman from FORT, maybe atoms

        self.dds_FORT.sw.on()
        self.dds_cooling_DP.sw.on()

        self.zotino0.set_dac(
            [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
            channels=self.coil_channels)

        self.ttl_SPCM_gate.on()
        delay(10*ms)

        if not self.disable_pulse_sequence:
            for j in range(1000):
                # self.laser_stabilizer.run()
                # self.dds_cooling_DP.sw.on()
                # self.dds_FORT.sw.on()

                self.core_dma.playback_handle(dma_handle)

                delay(1*ms)
        self.ttl_SPCM_gate.off()

