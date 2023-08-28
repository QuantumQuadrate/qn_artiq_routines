"""
Take an average with the sampler
"""

from artiq.experiment import *
import numpy as np

class SamplerAverage(EnvExperiment):

    def build(self):
        self.setattr_device("core")
        self.setattr_device("sampler0")
        self.setattr_argument("averages", NumberValue(4, type='int', ndecimals=0, scale=1, step=1))

    def prepare(self):
        self.data_array = np.zeros(8)
        self.measure_array = np.zeros(8)
        self.background_array = np.zeros(8)

    @rpc(flags={'async'})
    def print(self, average, data):
        print("average = ", average)
        print("data = ", data)

    @kernel
    def measure(self):
        for i in range(self.averages):
            print(i)
            delay(100*ms)
            self.sampler0.sample(self.data_array)
            self.measure_array += self.data_array
            delay(100*ms)
            self.print(self.measure_array, self.data_array)

        self.measure_array /= self.averages

    @kernel
    def background(self):
        for i in range(self.averages):
            print(i)
            delay(100*ms)
            self.sampler0.sample(self.data_array)
            self.background_array += self.data_array
            delay(100*ms)
            self.print(self.background_array, self.data_array)

        self.background_array /= self.averages

    @kernel
    def run(self):
        self.core.reset()

        # for i in range(self.averages):
        #     print(i)
        #     delay(100*ms)
        #     self.sampler0.sample(self.data_array)
        #     self.average_array += self.data_array
        #     delay(100*ms)
        #     self.print(self.average_array, self.data_array)
        #
        # self.average_array /= self.averages
        # self.print(self.average_array, self.data_array)

        # something like this gives issues in AOMPowerStabilizer
        self.measure()
        self.background()
        delay(100*ms)
        print(self.measure_array)
        delay(100*ms)
        print(self.measure_array - self.background_array)








