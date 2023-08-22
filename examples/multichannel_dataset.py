from artiq.experiment import *
import numpy as np
from numpy import pi
import time


class MultiChannelDataBroadcast(EnvExperiment):

    def build(self):
        self.setattr_argument("iterations", NumberValue(10, ndecimals=0, step=1))

    def run(self):

        # data_buffer is rolling: size is fixed, but we shift left and replace the last column
        y1 = [0.0] #np.zeros(1)
        y2 = [0.0] # np.zeros(1)
        pts = [10]

        self.set_dataset("channel1", y1, broadcast=True)
        self.set_dataset("channel2", y2, broadcast=True)
        self.set_dataset("pts", pts, broadcast=True) # how many pts to plot

        omega1 = 2*pi*np.random.rand()
        omega2 = 2*pi*np.random.rand()

        for i in range(self.iterations):
            a = np.random.rand()
            phi = 2*pi*np.random.rand()
            y1_next = a*np.sin(omega1*i + phi)
            a = np.random.rand()
            phi = 2 * pi * np.random.rand()
            y2_next = a * np.sin(omega2*i + phi)
            time.sleep(0.5)

            self.append_to_dataset("channel1", y1_next)
            self.append_to_dataset("channel2", y2_next)
            # self.append_to_dataset("pts", 0)