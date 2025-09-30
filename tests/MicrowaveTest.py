"""
Turns on the DDS channel which gets mixed with microwaves to drive the horn near 6.8 GHz.
The microwave switch is turned on to prevent leakage when microwave_dds_ON is false.
"""

from artiq.experiment import *

import sys, os
# get the current working directory
current_working_directory = os.getcwd()
cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")

from utilities.BaseExperiment import BaseExperiment


class MicrowaveTest(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.setattr_argument("microwave_dds_ON", BooleanValue(default=False))
        self.setattr_argument("microwave_dds_pulse_only", BooleanValue(default=False))
        self.setattr_argument("microwave_dds_pulse_time", NumberValue(10*ms, unit="ms"))

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

        # self.ttl7.pulse(5*ms)  # diagnostic trigger

        if not self.microwave_dds_pulse_only:
            if self.microwave_dds_ON:
                self.dds_microwaves.sw.on()
                self.ttl_microwave_switch.off()
            else:
                self.dds_microwaves.sw.off()
                self.ttl_microwave_switch.on()
        else:
            self.ttl_microwave_switch.off()
            self.dds_microwaves.sw.pulse(self.microwave_dds_pulse_time)
            self.ttl_microwave_switch.on()

        print("Experiment finished.")

