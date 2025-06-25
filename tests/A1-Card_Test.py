from artiq.experiment import *
from artiq.coredevice.ad9910 import (
    PHASE_MODE_ABSOLUTE,
    PHASE_MODE_CONTINUOUS,
    PHASE_MODE_TRACKING,
)

import numpy as np
# import pandas as pd
# import time


class Card_Tests(EnvExperiment):
# ### Testing TTLs:
#     def build(self):
#         self.setattr_device("core")
#         self.setattr_device("ttl13")
#         ### slef.setattr_device() is a built-in method in artiq that declares and configures a device such as core
#         ### and ttl4. These devices can be accessed like self.ttl4. //Akbar
#
#     @kernel
#     def run(self):
#         self.core.reset()
#
#         delay(1 * us)
#         self.ttl13.off()
#
#         for x in range(50):
#             self.ttl13.pulse(1*us)
#             delay(2*us)
#
#         print("code done!")



# ### Testing DDSs with phase control:
#     def build(self):
#        self.setattr_device("core")
#        self.setattr_device("urukul2_cpld")
#        self.setattr_device("urukul2_ch0")
#        self.setattr_device("urukul2_ch1")
#
#     @kernel
#     def run(self):
#         self.core.reset()
#         self.urukul2_ch0.cpld.init()
#         self.urukul2_ch0.init()
#         self.urukul2_ch0.set_att(float(0))
#
#         self.urukul2_ch1.cpld.init()
#         self.urukul2_ch1.init()
#         self.urukul2_ch1.set_att(float(0))
#
#         delay(1 * ms)
#
#         dBm = -10
#         self.urukul2_ch0.set_phase_mode(PHASE_MODE_TRACKING)
#         self.urukul2_ch1.set_phase_mode(PHASE_MODE_TRACKING)
#
#         self.urukul2_ch0.set(60.0 * MHz, amplitude=(2 * 50 * 10 ** (dBm / 10 - 3)) ** (1 / 2), phase = 0.0)
#         self.urukul2_ch1.set(60.0 * MHz, amplitude=(2 * 50 * 10 ** (dBm / 10 - 3)) ** (1 / 2), phase = 0.0)
#
#         self.urukul2_ch0.sw.on()
#         self.urukul2_ch1.sw.on()
#
#         delay(5000 * ms)
#         self.urukul2_ch1.set(60.0 * MHz, amplitude=(2 * 50 * 10 ** (dBm / 10 - 3)) ** (1 / 2), phase=0.5)
#
#         # self.urukul2_ch0.sw.off()
#         # self.urukul2_ch1.sw.off()
#
#         print("code done!")


### Testing DDSs:
    def build(self):
       self.setattr_device("core")
       # print("aaaa")
       self.setattr_device("urukul2_cpld")
       self.setattr_device("urukul2_ch3")


    @kernel
    def run(self):
        self.core.reset()
        self.urukul2_ch3.cpld.init()
        self.urukul2_ch3.init()
        self.urukul2_ch3.set_att(float(0))

        # self.urukul0_cpld.set_profile(0)

        delay(10 * ms)

        dBm = -12
        self.urukul2_ch3.set(334.682 * MHz, amplitude=(2 * 50 * 10 ** (dBm / 10 - 3)) ** (1 / 2)) #0.08)
        self.urukul2_ch3.sw.on()

        delay(2000 * ms)
        self.urukul2_ch3.sw.off()

        print("code done!")


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
 #       self.zotino0.write_dac(0, 0.0)
 #       self.zotino0.write_dac(1, 0.0)
 #       self.zotino0.write_dac(2, 0.0)
 #       self.zotino0.write_dac(3, 0.0)
 #       self.zotino0.write_dac(4, 0.0)
 #       self.zotino0.write_dac(5, 0.0)
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






#### reading Sampler - with averaging:
    # def build(self):
    #     self.setattr_device("core")
    #     self.setattr_device("sampler0")
    #
    # def prepare(self):
    #     n_channels = 8
    #     self.smp = np.zeros(n_channels, dtype=float)
    #     self.avg = np.zeros(n_channels, dtype=float)
    #
    # @kernel
    # def run(self):
    #     self.core.reset()
    #     self.core.break_realtime()
    #     self.sampler0.init()
    #
    #     for i in range(8):  # loops for each sampler channel
    #        self.sampler0.set_gain_mu(i, 0)  # sets each channel's gain to 0db
    #        delay(100 * us)  # 100us delay
    #
    #     iters = 100
    #
    #     for i in range(iters):
    #         self.sampler0.sample(self.smp)  # runs sampler and saves to list
    #         self.avg += self.smp
    #         delay(0.1*ms)
    #     self.avg /= iters
    #
    #     print("sampler0 channel-wise average")
    #     print(self.avg)





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






