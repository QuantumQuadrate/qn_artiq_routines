"""
Turns on the DDS channel which gets mixed with microwaves to drive the horn near 6.8 GHz.
The microwave switch is turned on to prevent leakage when microwave_dds_ON is false.
"""

from artiq.experiment import *

import sys, os
sys.path.append('C:\\Users\\jakeuribe\\artiq-master\\repository\\qn_artiq_routines\\')
sys.path.append('C:\\Users\\jakeuribe\\artiq-master\\repository\\')


from utilities.BaseExperiment import BaseExperiment


class MicrowaveTest(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.setattr_argument("microwave_dds_ON", BooleanValue(default=False))

        self.base.set_datasets_from_gui_args()
        print("build - done")

    def prepare(self):
        """
        performs initial calculations and sets parameter values before
        running the experiment. also sets data filename for now.

        any conversions from human-readable units to machine units (mu) are done here
        """
        self.base.prepare()

        print("prepare - done")

    @kernel
    def run(self):
        self.base.initialize_hardware()
        self.expt()

        if self.microwave_dds_ON:
            self.dds_microwaves.sw.on()
            self.ttl_microwave_switch.off()
        else:
            self.dds_microwaves.sw.off()
            self.ttl_microwave_switch.on()

        print("Experiment finished.")

    @kernel
    def expt(self):
        """
        The experiment loop.

        :return:
        """
