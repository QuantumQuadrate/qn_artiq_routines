from artiq.experiment import *
from artiq.coredevice.urukul import *


class Urukul0Test(EnvExperiment):

    def build(self):
        self.setattr_device('core')
        self.setattr_device("urukul0_cpld")
        self.setattr_device("urukul0_ch0")
        self.setattr_device("urukul0_ch1")
        self.setattr_device("urukul0_ch2")
        self.setattr_device("urukul0_ch3")

    def prepare(self):
        pass

    @kernel
    def run(self):
        self.core.reset()
        self.urukul0_cpld.init()  # proto_rev mismatch on Node 1
        self.urukul0_ch0.init()
        self.urukul0_ch1.init()
        self.urukul0_ch2.init()
        self.urukul0_ch3.init()
        self.urukul0_ch0.set_att(float(0))
        self.urukul0_ch1.set_att(float(0))
        self.urukul0_ch2.set_att(float(0))
        self.urukul0_ch3.set_att(float(0))

        self.urukul0_ch0.sw.on()
        delay(0.5 * s)
        self.urukul0_ch0.sw.off()
        self.urukul0_ch1.sw.on()
        delay(0.5 * s)
        self.urukul0_ch1.sw.off()
        self.urukul0_ch2.sw.on()
        delay(0.5 * s)
        self.urukul0_ch2.sw.off()
        self.urukul0_ch3.sw.on()
        delay(0.5 * s)
        self.urukul0_ch3.sw.off()
