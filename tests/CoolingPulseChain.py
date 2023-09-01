"""
Chop the cooling double pass AOM. Can be used for measuring detector response.
"""

from artiq.experiment import *

from utilities.BaseExperiment import BaseExperiment


class CoolingPulseChain(EnvExperiment):
    def build(self):
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.setattr_argument("n_cycles", NumberValue(1000, type='int', scale=1, ndecimals=0, step=1))
        self.setattr_argument("period", NumberValue(1.0, unit='ms', ndecimals=2))

    def prepare(self):
        self.base.prepare()
        self.half_period = self.period/2

    @kernel
    def run(self):
        self.base.initialize_hardware()

        for i in range(self.n_cycles):
            self.dds_cooling_DP.sw.on()
            delay(self.half_period)
            self.dds_cooling_DP.sw.off()
            delay(self.half_period)

        print("Cooling pulse chain done!")