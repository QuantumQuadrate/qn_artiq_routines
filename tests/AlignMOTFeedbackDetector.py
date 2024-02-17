"""
Repeatedly loop over the fiber AOMs 1-4, turning only one on at a time.
This allows us to position a detector to try to equalize the signal from each.
"""

from artiq.experiment import *

import sys, os
sys.path.append('C:\\Users\\jakeuribe\\artiq-master\\repository\\qn_artiq_routines\\')
sys.path.append('C:\\Users\\jakeuribe\\artiq-master\\repository\\')


from utilities.BaseExperiment import BaseExperiment

class AlignMOTFeedbackDetector(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.setattr_argument("run_time_minutes", NumberValue(1))
        self.setattr_argument("dt_AOM_on", NumberValue(10 * ms, unit='ms'))

        self.base.set_datasets_from_gui_args()

        print("build - done")

    def prepare(self):
        """
        performs initial calculations and sets parameter values before
        running the experiment. also sets data filename for now.

        any conversions from human-readable units to machine units (mu) are done here
        """

        self.base.prepare()
        self.n_steps = int(60*self.run_time_minutes/(4*self.dt_AOM_on)+0.5)

        print("prepare - done")

    @kernel
    def run(self):
        self.base.initialize_hardware()

        self.dds_cooling_DP.sw.on()

        self.dds_AOM_A2.sw.on()
        self.dds_AOM_A3.sw.on()
        self.dds_AOM_A1.sw.on()
        self.dds_AOM_A4.sw.on()

        # wait for AOMs to thermalize
        delay(3000 * ms)

        print("ready!")

        for i in range(self.n_steps):

            self.dds_AOM_A1.sw.on()
            self.dds_AOM_A2.sw.off()
            self.dds_AOM_A3.sw.off()
            self.dds_AOM_A4.sw.off()
            delay(self.dt_AOM_on)

            self.dds_AOM_A1.sw.off()
            self.dds_AOM_A2.sw.on()
            self.dds_AOM_A3.sw.off()
            self.dds_AOM_A4.sw.off()
            delay(self.dt_AOM_on)

            self.dds_AOM_A1.sw.off()
            self.dds_AOM_A2.sw.off()
            self.dds_AOM_A3.sw.on()
            self.dds_AOM_A4.sw.off()
            delay(self.dt_AOM_on)

            self.dds_AOM_A1.sw.off()
            self.dds_AOM_A2.sw.off()
            self.dds_AOM_A3.sw.off()
            self.dds_AOM_A4.sw.on()
            delay(self.dt_AOM_on)

            self.dds_AOM_A1.sw.off()
            self.dds_AOM_A2.sw.off()
            self.dds_AOM_A3.sw.off()
            self.dds_AOM_A4.sw.off()
            delay(self.dt_AOM_on)

        print("Experiment finished.")
