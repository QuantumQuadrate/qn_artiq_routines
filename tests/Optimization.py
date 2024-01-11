from artiq.experiment import *
import numpy as np
import scipy as sp

class Optimization(EnvExperiment):

    def build(self):
        self.setattr_device("core")
        self.setattr_device("led0")

    def prepare(self):

        x_pts = np.linspace(0,100,100) - 11
        self.xmean = 11
        self.ymean = 20
        domain = 100 # the constrained space
        self.xmin = - domain/2 - self.xmean
        self.xmax = domain/2 + self.xmean
        self.ymin = - domain / 2 - self.ymean
        self.ymax = domain / 2 + self.ymean

    def landscape(self,arr: TArray(TFloat,1)) -> TFloat:
        """some surface that we want to minimize"""

        return -1*np.exp(-((arr[0]-self.xmean)/20)**2)*np.exp(-((arr[1]-self.ymean)/10)**2)

    @kernel
    def experiment(self, arr: TArray(TFloat,1)) -> TFloat:
        self.core.reset()
        self.led0.pulse(10*us) # just because
        return self.landscape(arr) # we are trying to maximize landscape so negate it for the minimizer

    def run(self):

        result = sp.optimize.minimize(fun=self.experiment,
                                      x0=[11.0,10],
                                      bounds=[[self.xmin,self.xmax],[self.ymin,self.ymax]],
                                      method='Nelder-Mead')
        print("optimization successful?",result.success)
        print("solution:",result.x)

        print("experiment finished")