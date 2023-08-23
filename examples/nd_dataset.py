from artiq.experiment import *
import numpy as np
from numpy import pi
import time


class NDDataBroadcast(EnvExperiment):

    def build(self):
        self.setattr_argument("iterations", NumberValue(10, ndecimals=0, step=1))
        self.setattr_argument("channels", NumberValue(5, ndecimals=0, step=1))

    def run(self):

        # data_buffer is rolling: size is fixed, but we shift left and replace the last column
        data_buffer = np.full((self.channels, self.iterations), 0.0)

        self.set_dataset("nd_data", data_buffer, broadcast=True)

        a = [np.random.rand() for i in range(self.channels)]
        phi = [2*pi*np.random.rand() for i in range(self.channels)]
        omega = [2*pi*np.random.rand() for i in range(self.channels)]

        latest_data = np.array([0.0]*self.channels)

        for i in range(self.iterations):
            for ch in range(self.channels):
                latest_data[ch] = (0.5+0.5*np.random.rand())*a[ch]*np.sin(omega[ch]*i + phi[ch])
            time.sleep(0.5)
            print(latest_data)

            data_buffer = np.roll(data_buffer, -1, 1) # shift columns to the left
            data_buffer[:,-1] = latest_data
            print(data_buffer)

            self.mutate_dataset("nd_data", (0,self.channels), data_buffer)