from artiq.experiment import *
from artiq.coredevice.urukul import *


class GetProtoRev(EnvExperiment):

    def build(self):
        self.setattr_device('core')
        self.setattr_device("core")
        self.setattr_device("urukul0_cpld")
        self.setattr_device("urukul1_cpld")
        self.setattr_device("urukul2_cpld")

    def prepare(self):
        pass

    @kernel
    def run(self):
        self.core.reset()
        # self.urukul0_cpld.init() # proto_rev mismatch on Node 1

        proto_rev = urukul_sta_proto_rev(self.urukul0_cpld.sta_read())
        print("cpld0 proto_rev: ", proto_rev)
        print("expected:", STA_PROTO_REV_MATCH)

        # delay(0.5*s)
        # proto_rev = urukul_sta_proto_rev(self.urukul1_cpld.sta_read())
        # print("cpld1 proto_rev: ", proto_rev)
        # print("expected:", STA_PROTO_REV_MATCH)
        #
        # delay(0.5*s)
        # proto_rev = urukul_sta_proto_rev(self.urukul2_cpld.sta_read())
        # print("cpld2 proto_rev: ", proto_rev)
        # print("expected:", STA_PROTO_REV_MATCH)
