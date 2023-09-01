"""
Pulse an AOM to measure the rise time
"""

from artiq.experiment import *
import math

import sys, os
sys.path.append('C:\\Networking Experiment\\artiq codes\\artiq-master\\repository\\qn_artiq_routines\\')

from utilities.BaseExperiment import BaseExperiment

class FiberAOMRiseTimeTest(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.setattr_argument("run_time_minutes", NumberValue(1))
        self.setattr_argument("dt_AOM_on", NumberValue(1000 * ms, unit='ms'))
        self.setattr_argument("which_fiber_AOM", NumberValue(1, type='int', ndecimals=0, scale=1, step=1))
        self.setattr_argument("leave_on_but_change_frequency", BooleanValue(False))
        self.setattr_argument("default_frequency", NumberValue(78.5 * MHz, unit='MHz'))
        self.setattr_argument("default_power", NumberValue(-8, unit="dBm",scale=1, ndecimals=1))
        self.setattr_argument("switch_frequency", NumberValue(70 * MHz, unit='MHz'))

        print("build - done")

    def prepare(self):
        """
        performs initial calculations and sets parameter values before
        running the experiment. also sets data filename for now.

        any conversions from human-readable units to machine units (mu) are done here
        """

        self.base.prepare()
        self.n_steps = int(60*self.run_time_minutes/(4*self.dt_AOM_on)+0.5)
        AOM_dict = {1: self.dds_AOM_A1,
                    2: self.dds_AOM_A2,
                    3: self.dds_AOM_A3,
                    4: self.dds_AOM_A4,
                    5: self.dds_AOM_A5,
                    6: self.dds_AOM_A6}
        self.AOM = AOM_dict[self.which_fiber_AOM]
        self.default_amplitude = math.sqrt(2*50*10**(self.default_power/10-3))

        print("prepare - done")

    @kernel
    def run(self):
        self.base.initialize_hardware()

        delay(10*ms)

        for i in range(self.n_steps):
            if self.leave_on_but_change_frequency:
                self.AOM.set(frequency=self.default_frequency, amplitude=self.default_amplitude)
                delay(self.dt_AOM_on)
                self.AOM.set(frequency=self.switch_frequency, amplitude=self.default_amplitude)
                delay(self.dt_AOM_on)

            else:
                self.AOM.sw.on()
                delay(self.dt_AOM_on)
                self.AOM.sw.off()
                delay(self.dt_AOM_on)

        self.AOM.set(frequency=self.default_frequency, amplitude=self.default_amplitude)
        self.AOM.sw.off()

        print("Experiment finished.")
