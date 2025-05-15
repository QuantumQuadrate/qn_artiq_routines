"""
Testing how dds profiles work, how I can set them and call them. Works well. Can change the profile as quickly as 10ns without issue.
However, the first profile stay up for about 200ns.

Akbar 2025-05-06

"""
from artiq.experiment import *

class DDS_profile_test(EnvExperiment):
    def build(self):
        self.setattr_device("core")
        self.setattr_device("urukul2_cpld")
        self.setattr_device("urukul2_ch0")

    @kernel
    def run(self):
        self.core.reset()
        self.urukul2_cpld.init()
        self.urukul2_ch0.init()
        self.urukul2_ch0.set_att(0.0)

        ### Load profile 0
        self.urukul2_ch0.set(100.0*MHz, amplitude=1.0, profile=0)

        ### Load profile 1
        self.urukul2_ch0.set(150.0*MHz, amplitude=0.5, profile=1)

        delay(10*ms)

        ### Select profile 0
        self.urukul2_cpld.set_profile(0)
        self.urukul2_cpld.io_update.pulse_mu(8)
        self.urukul2_ch0.sw.on()
        delay(10*ns)

        ### Switch to profile 1
        self.urukul2_cpld.set_profile(1)
        self.urukul2_cpld.io_update.pulse_mu(8)
        delay(5*ms)

        self.urukul2_ch0.sw.off()



