"""
Take an average with the sampler
"""

from artiq.experiment import *
import numpy as np

class SamplerAverage(EnvExperiment):

    def build(self):
        self.setattr_device("core")
        self.setattr_device("sampler0")
        self.setattr_argument("averages", NumberValue(4, type='int', precision=0, scale=1, step=1))

    def prepare(self):
        self.data_array = np.zeros(8)
        self.measure_array = np.zeros(8)
        self.background_array = np.zeros(8)

    @rpc(flags={'async'})
    def print(self, x):
        print(x)

    @kernel
    def measure(self):
        self.measure_array = np.full(8, 0.0)
        for i in range(self.averages):
            self.sampler0.sample(self.data_array)
            self.measure_array += self.data_array
            delay(1*ms)
            # self.print(self.measure_array, self.data_array)

        self.measure_array /= self.averages
        self.print("self.measure_array")
        self.print(self.measure_array)

    @kernel
    def background(self):
        self.background_array = np.full(8, 0.0)
        for i in range(self.averages):
            self.sampler0.sample(self.data_array)
            self.background_array += self.data_array
            delay(1*ms)

        self.background_array /= self.averages
        self.print("self.background_array")
        self.print(self.background_array)

    @kernel
    def run(self):
        self.core.reset()

        delay(1*ms)
        self.background()
        delay(1*ms)
        self.measure()
        delay(1*ms)
        self.background()







