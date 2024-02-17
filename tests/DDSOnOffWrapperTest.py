"""
The microwave test code but with a wrapper for the dds's switch on,off methods

adding attributes to the urukul object in prepare doesn't succeed, and I doubt i can do that in run
"""

from artiq.experiment import *

import sys
sys.path.append('C:\\Users\\jakeuribe\\artiq-master\\repository\\qn_artiq_routines\\')
sys.path.append('C:\\Users\\jakeuribe\\artiq-master\\repository\\')


from utilities.BaseExperiment import BaseExperiment


class DDSOnOffWrapperTest(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        # self.setattr_argument("microwave_dds_ON", BooleanValue(default=False))

        self.base.set_datasets_from_gui_args()
        print("build - done")

    def prepare(self):
        """
        performs initial calculations and sets parameter values before
        running the experiment. also sets data filename for now.

        any conversions from human-readable units to machine units (mu) are done here
        """
        self.base.prepare()

        def on_wrapper(sw):
            @kernel
            def wrapper():
                sw.on_original()
                sw.is_on = True
            return wrapper

        def off_wrapper(sw):
            @kernel
            def wrapper():
                sw.off_original()
                sw.is_on = False
            return wrapper

        # self.dds_microwaves.sw.on_original = self.dds_microwaves.sw.on
        # self.dds_microwaves.sw.off_original = self.dds_microwaves.sw.off
        # self.dds_microwaves.sw.is_on = False
        # self.dds_microwaves.sw.on = on_wrapper(self.dds_microwaves.sw)
        # self.dds_microwaves.sw.off = off_wrapper(self.dds_microwaves.sw)
        setattr(self.dds_microwaves.sw, "on_original", self.dds_microwaves.sw.on)
        setattr(self.dds_microwaves.sw, "off_original", self.dds_microwaves.sw.off)
        setattr(self.dds_microwaves.sw, "is_on", False)
        setattr(self.dds_microwaves.sw, "on", on_wrapper(self.dds_microwaves.sw))
        setattr(self.dds_microwaves.sw, "off", off_wrapper(self.dds_microwaves.sw))

        print("prepare - done")

    @kernel
    def run(self):
        self.base.initialize_hardware()
        self.expt()

        self.print_async("turning microwaves dds channel on")
        self.dds_microwaves.sw.on()
        # self.ttl_microwave_switch.off()
        self.print_async("dds on =",self.dds_microwaves.sw.is_on)

        delay(1*s)

        self.print_async("turning microwaves dds channel off")
        self.dds_microwaves.sw.off()
        # self.ttl_microwave_switch.on()
        self.print_async("dds on =", self.dds_microwaves.sw.is_on)

        print("Experiment finished.")

    @kernel
    def expt(self):
        """
        The experiment loop.

        :return:
        """
