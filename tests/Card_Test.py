from artiq.experiment import *
import numpy as np
# import pandas as pd
# import time


class Card_Tests(EnvExperiment):
### Testing TTLs:
    def build(self):
        self.setattr_device("core")
        self.setattr_device("ttl4")
        ### slef.setattr_device() is a built-in method in artiq that declares and configures a device such as core
        ### and ttl4. These devices can be accessed like self.ttl4. //Akbar

    @kernel
    def run(self):
        self.core.reset()
        self.ttl4.output()

        delay(1 * us)

        # self.ttl4.on()
        self.ttl4.off()

        # for x in range(50):
        #     self.ttl4.pulse(10*ms)
        #     delay(10*ms)

        print("code done!")


# #### Testing DDSs:
#     def build(self):
#        self.setattr_device("core")
#        self.setattr_device("urukul0_cpld")
#        self.setattr_device("urukul0_ch0")
#        self.setattr_device("urukul0_ch3")
#
#     @kernel
#     def run(self):
#        self.core.reset()
#        self.urukul0_ch3.cpld.init()
#        self.urukul0_ch3.init()
#        self.urukul0_ch3.set_att(float(0))
#
#        # self.urukul0_cpld.set_profile(0)
#
#        delay(10 * ms)
#
#        self.urukul0_ch3.set(150.5 * MHz, amplitude=0.9)
#        self.urukul0_ch3.sw.on()
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
#        self.setattr_device("sampler0")
#
#     @kernel
#     def run(self):
#        self.core.reset()
#        self.core.break_realtime()
#        self.sampler0.init()
#
#        for i in range(8):  # loops for each sampler channel
#            self.sampler0.set_gain_mu(i, 0)  # sets each channel's gain to 0db
#            delay(100 * us)  # 100us delay
#
#        n_channels = 8
#
#        delay(10 * ms)
#
#        smp = [0.0] * n_channels
#
#        self.sampler0.sample(smp)  # runs sampler and saves to list
#
#        for i in range(len(smp)):  # loops over list of samples
#            print("ch",i,":",smp[i])
#            delay(100*ms)
# #
# #        delay(10 * ms)
# #
# #        print("code done!")






# #### reading Sampler - with averaging:
#     def build(self):
#         self.setattr_device("core")
#         self.setattr_device("sampler0")
#
#     def prepare(self):
#         n_channels = 8
#         self.smp = np.zeros(n_channels, dtype=float)
#         self.avg = np.zeros(n_channels, dtype=float)
#
#     @kernel
#     def run(self):
#         self.core.reset()
#         self.core.break_realtime()
#         self.sampler0.init()
#
#         for i in range(8):  # loops for each sampler channel
#            self.sampler0.set_gain_mu(i, 0)  # sets each channel's gain to 0db
#            delay(100 * us)  # 100us delay
#
#         iters = 100
#
#         for i in range(iters):
#             self.sampler0.sample(self.smp)  # runs sampler and saves to list
#             self.avg += self.smp
#             delay(0.1*ms)
#         self.avg /= iters
#
#         print("sampler0 channel-wise average")
#         print(self.avg)





#### Sending a series of triggering signals to the Thorcam. Since the TTL signal is too fast,
#### we are using zotino
    # def build(self):
    #     self.setattr_device("core")
    #     self.setattr_device("zotino0")
    #
    # @kernel
    # def run(self):
    #     self.core.reset()
    #     self.core.break_realtime()
    #     self.zotino0.init()
    #     delay(10 * ms)
    #
    #     self.zotino0.write_dac(6, 4.0)
    #     self.zotino0.load()
    #     delay(5 * ms)
    #     self.zotino0.write_dac(6, 0.0)
    #     self.zotino0.load()
    #
    #
    #     delay(10 * ms)
    #
    #     print("code done!")






