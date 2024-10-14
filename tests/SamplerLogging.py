"""
This code measures analog signals with the Samplers, save them in a dataset to plot live and for logging.
"""
from artiq.experiment import *
import numpy as np


class SamplerLogging(EnvExperiment):

    def build(self):
        self.setattr_device("core")
        self.setattr_device("sampler0")
        self.setattr_device("sampler1")
        self.setattr_device("sampler2")

        self.setattr_argument("n_average", NumberValue(5, type='int', ndecimals=0, scale=1, step=1))  # averaging over n in each measurement
        self.setattr_argument("n_measure", NumberValue(10, type='int', ndecimals=0, scale=1, step=1))  # number of measurements
        self.setattr_argument("t_step", NumberValue(10, type='int', unit='ms', ndecimals=0, scale=1, step=1))  # delay between measurements
        self.setattr_argument("t_step_in_average", NumberValue(1, type='int', unit='ms', ndecimals=0, scale=1, step=1))  # delay between measurements

    def prepare(self):
        self.n_channels = 8
        self.smp = np.zeros(self.n_channels, dtype=float)
        self.avg = np.zeros(self.n_channels, dtype=float)
        self.smpList = [self.avg]


    @kernel
    def run(self):
        self.core.reset()
        self.core.break_realtime()
        self.sampler0.init()
        self.sampler1.init()
        self.sampler2.init()

        self.sampler0.sample(self.smp)  # runs sampler and saves to list
        delay(1 * ms)
        self.set_dataset("SamplerValues", self.smpList, broadcast=True, persist=True)
        self.set_dataset("t_step", self.t_step, broadcast=True, persist=True)
        delay(1 * ms)


        for i in range(self.n_measure):
            delay(self.t_step * ms)
            for j in range(self.n_average):
                self.sampler0.sample(self.smp)  # runs sampler and saves to list
                self.avg += self.smp
                delay(self.t_step_in_average * ms)
            self.avg /= self.n_average

            self.append_to_dataset("SamplerValues", self.avg)






