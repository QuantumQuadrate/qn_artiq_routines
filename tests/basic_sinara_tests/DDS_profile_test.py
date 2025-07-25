"""
Testing how dds profiles work, how I can set them and call them. Works well. Can change the profile as quickly as 10ns without issue.
When using profile, the delay from the time we switch the profile until it is applied and the amplitude, for example, changes, is about
200ns (observed from oscilloscope). However, when not using profile (changing the amplitude by dds.set(...) as shown below) this delay is
about 2us. Other than that, I did not notice any difference between the two methods and the waveform on the scope is continuous in
both cases.

Akbar 2025-05-06

"""
from artiq.experiment import *

class DDS_profile_test(EnvExperiment):
    def build(self):
        self.setattr_device("core")
        self.setattr_device("ttl7")
        self.dds1 = self.get_device("urukul2_ch0")
        self.dds2 = self.get_device("urukul2_ch1")

    @kernel
    def run(self):
        self.core.reset()

        self.dds1.cpld.init()
        self.dds1.init()
        self.dds1.set_att(0.0)

        self.dds2.cpld.init()
        self.dds2.init()
        self.dds2.set_att(0.0)

        self.using_profile()
        # self.without_profile()


    @kernel
    def using_profile(self):
        ### Load profile 0
        self.dds1.set(10.0 * MHz, amplitude=0.2, profile=0)
        self.dds2.set(10.0 * MHz, amplitude=0.5, profile=0)

        ### Load profile 1
        self.dds1.set(10.0 * MHz, amplitude=0.4, profile=1)
        self.dds2.set(10.0 * MHz, amplitude=0.5, profile=1)
        ### Even if we don't want to change dd2, since it is on the same card as dds1, we need to define profile=1 for that too. Otherwise,
        ### when switching to profile 1, dds2 will turn off.

        delay(10 * ms)

        self.ttl7.on()
        self.dds2.sw.on()

        ### Select profile 0
        self.dds1.cpld.set_profile(0)
        self.dds1.cpld.io_update.pulse_mu(8)
        self.dds1.sw.on()
        delay(2000 * ns)

        ### Switch to profile 1
        self.dds1.cpld.set_profile(1)
        self.dds1.cpld.io_update.pulse_mu(8)
        delay(5 * us)

        self.dds1.sw.off()
        self.dds2.sw.off()

        self.ttl7.off()

    @kernel
    def without_profile(self):
        delay(1 * ms)

        self.dds1.set(10.0 * MHz, amplitude=0.2)
        self.dds2.set(10.0 * MHz, amplitude=0.5)

        self.ttl7.on()
        self.dds1.sw.on()
        self.dds2.sw.on()

        delay(2000 * ns)

        self.dds1.set(10.0 * MHz, amplitude=0.4)
        self.dds2.set(10.0 * MHz, amplitude=0.5)

        delay(5 * us)

        self.dds1.sw.off()
        self.dds2.sw.off()

        self.ttl7.off()






