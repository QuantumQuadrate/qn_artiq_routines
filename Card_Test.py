from artiq.experiment import *
# import pandas as pd
# import time

class Card_Tests(EnvExperiment):

# ### Testing TTLs:
#     def build(self):
#         self.setattr_device("core")
#         self.setattr_device("ttl4")
#
#     @kernel
#     def run(self):
#         self.core.reset()
#         self.ttl4.output()
#
#         delay(1 * us)
#
#         for x in range(500):
#             self.ttl4.pulse(10*ms)
#             delay(10*ms)
#
#         print("code done!")


# #### Testing DDSs:
#     def build(self):
#        self.setattr_device("core")
#        self.setattr_device("urukul0_cpld")
#        self.setattr_device("urukul0_ch0")
#        self.setattr_device("urukul0_ch1")
#
#     @kernel
#     def run(self):
#        self.core.reset()
#        self.urukul0_ch0.cpld.init()
#        self.urukul0_ch0.init()
#        self.urukul0_ch0.set_att(float(0))
#
#        # self.urukul0_cpld.set_profile(0)
#
#        delay(10 * ms)
#
#        self.urukul0_ch0.set(70.0 * MHz, amplitude=0.9)
#        self.urukul0_ch0.sw.on()
#
#        # delay(2000 * ms)
#        # self.urukul0_ch0.sw.off()
#
#        print("code done!")


 # #### Testing Zotino:
 #    def build(self):
 #       self.setattr_device("core")
 #       self.setattr_device("zotino0")
 #
 #    @kernel
 #    def run(self):
 #       self.core.reset()
 #       self.core.break_realtime()
 #       self.zotino0.init()
 #       delay(10 * ms)
 #
 #       self.zotino0.write_dac(0,1.0)
 #       self.zotino0.load()
 #
 #       # self.zotino0.set_dac([1.0], [6])
 #       # self.zotino0.set_dac([2,2,2,2])
 #
 #       delay(10 * ms)
 #
 #       print("code done!")



# #### Testing Sampler:
#     def build(self):
#        self.setattr_device("core")
#        self.setattr_device("sampler2")
#
#     @kernel
#     def run(self):
#        self.core.reset()
#        self.core.break_realtime()
#        self.sampler2.init()
#
#        for i in range(8):  # loops for each sampler channel
#            self.sampler2.set_gain_mu(i, 0)  # sets each channel's gain to 0db
#            delay(100 * us)  # 100us delay
#
#        n_channels = 8
#
#        delay(10 * ms)
#
#        smp = [0.0] * n_channels
#
#        self.sampler2.sample(smp)  # runs sampler and saves to list
#
#        for i in range(len(smp)):  # loops over list of samples
#            print("ch",i,":",smp[i])
#            delay(100*ms)
# #
# #        delay(10 * ms)
# #
# #        print("code done!")




#### Sending a series of triggering signals to the Thorcam. Since the TTL signal is too fast,
#### we are using zotino
    def build(self):
        self.setattr_device("core")
        self.setattr_device("zotino0")

    @kernel
    def run(self):
        self.core.reset()
        self.core.break_realtime()
        self.zotino0.init()
        delay(10 * ms)

        self.zotino0.write_dac(6, 4.0)
        self.zotino0.load()
        delay(5 * ms)
        self.zotino0.write_dac(6, 0.0)
        self.zotino0.load()


        delay(10 * ms)

        print("code done!")





