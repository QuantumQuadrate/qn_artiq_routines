"""
When a photon gets detected by the SPCM, and we want to trigger the microwave pulse for atom mapping and rotation,
is the delay between the SPCM TTL output and the microwave signal deterministic?

I wrote this simple code to test this. Use a T to split the SPCM output one goes to TTL count one goes to a scope
for triggering. The TTL out for microwave switch also goes to the same scope. I confirmed that the delay between
these two signals is fixed with about 1ns jitter. With self.t_start_MW_mapping_mu = 1000 I was getting Underflow
error ocassinally, but it was good enough to run this experiment. The SPCM is just detecting noises.

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


class Card_Tests(EnvExperiment):
### Testing TTLs:
    def build(self):
        self.base = BaseExperiment(experiment=self)
        self.base.build()

    def prepare(self):
        self.base.prepare()

    @kernel
    def run(self):
        self.core.reset()

        delay(1 * us)
        self.ttl_microwave_switch.off()
        delay(1 * ms)
        self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves))
        delay(10*us)
        self.dds_microwaves.sw.on()

        self.t_start_MW_mapping_mu = 1000

        for x in range(10000):
            t_end_SPCM0 = self.ttl_SPCM0.gate_rising(1*ms)
            SPCM0_click_time = self.ttl_SPCM0.timestamp_mu(t_end_SPCM0)

            if SPCM0_click_time > 0:
                at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu)

                self.ttl_microwave_switch.off()
                delay(10*us)
                self.ttl_microwave_switch.on()

            delay(1*ms)


        print("code done!")
