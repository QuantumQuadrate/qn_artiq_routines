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
        self.setattr_argument("t_step", NumberValue(10, type='float', unit='ms', ndecimals=2, scale=1, step=1))  # delay between measurements
        self.setattr_argument("t_step_in_average", NumberValue(1, type='float', unit='ms', ndecimals=2, scale=1, step=1))  # delay between measurements

    def prepare(self):
        self.n_channels = 8

        self.smp0 = np.zeros(self.n_channels, dtype=float)
        self.smp1 = np.zeros(self.n_channels, dtype=float)
        self.smp2 = np.zeros(self.n_channels, dtype=float)

        self.avg0 = np.zeros(self.n_channels, dtype=float)
        self.avg1 = np.zeros(self.n_channels, dtype=float)
        self.avg2 = np.zeros(self.n_channels, dtype=float)

        self.smp0List = [self.avg0]
        self.smp1List = [self.avg1]
        self.smp2List = [self.avg2]


    @rpc(flags={"async"})
    def print_async(*x):
        """print asynchronously so we don't block the RTIO counter.
        useful for debugging"""
        print(*x)


    @kernel
    def run(self):
        self.core.reset()
        self.core.break_realtime()
        self.sampler0.init()
        self.sampler1.init()
        self.sampler2.init()

        self.set_dataset("Sampler0Values", self.smp0List, broadcast=True, persist=False)
        self.set_dataset("Sampler1Values", self.smp1List, broadcast=True, persist=False)
        self.set_dataset("Sampler2Values", self.smp2List, broadcast=True, persist=False)
        self.set_dataset("t_step", self.t_step, broadcast=True, persist=False)

        # for i in range(8):
        #     self.sampler0.set_gain_mu(i, 10)
        #     delay(100 * us)

        progress_steps = 10
        step = max(1, self.n_measure // progress_steps)

        for i in range(self.n_measure):
            dummy0 = np.full(8, 0.0)
            dummy1 = np.full(8, 0.0)
            dummy2 = np.full(8, 0.0)
            delay(self.t_step * ms)
            for j in range(self.n_average):
                self.sampler0.sample(self.smp0)  # runs sampler and saves to list
                delay(0.1 * ms)
                self.sampler1.sample(self.smp1)  # runs sampler and saves to list
                delay(0.1 * ms)
                self.sampler2.sample(self.smp2)  # runs sampler and saves to list
                delay(0.1 * ms)
                dummy0 += self.smp0
                dummy1 += self.smp1
                dummy2 += self.smp2
                delay(self.t_step_in_average * ms)

            dummy0 /= self.n_average
            dummy1 /= self.n_average
            dummy2 /= self.n_average

            self.avg0 = dummy0
            self.avg1 = dummy1
            self.avg2 = dummy2

            self.append_to_dataset("Sampler0Values", self.avg0)
            self.append_to_dataset("Sampler1Values", self.avg1)
            self.append_to_dataset("Sampler2Values", self.avg2)

            delay(1*ms)
            if i % step == 0 or i == self.n_measure - 1:
                self.print_async("Progress: ", 100 * i / self.n_measure)





