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
#         self.setattr_device("ttl15")
#
#     @kernel
#     def run(self):
#         self.core.reset()
#
#         delay(1 * us)
#         self.ttl15.off()
#
#         for x in range(1000):
#             self.ttl15.pulse(1*ms)
#             delay(100*ms)
#
#         # self.ttl15.off()
#
#         print("code done!")




### Testing DDSs with phase control:
### Without phase control, the phase between different dds channels is still fixed over time, even if the channels
### are on different cards. Adding a delay (turning a channel on with some delay after the other channel) does not
### affect the phase relation between different channels.

    # def build(self):
    #    self.setattr_device("core")
    #    self.setattr_device("urukul2_cpld")
    #    self.setattr_device("urukul2_ch0")
    #    self.setattr_device("urukul2_ch1")
    #
    # @kernel
    # def run(self):
    #     self.core.reset()
    #     self.urukul2_cpld.init()
    #     self.urukul2_ch0.init()
    #     self.urukul2_ch0.set_att(float(0))
    #
    #     self.urukul2_ch1.init()
    #     self.urukul2_ch1.set_att(float(0))
    #
    #     delay(1 * ms)
    #
    #     dBm = -5
    #     # self.urukul2_ch0.set_phase_mode(PHASE_MODE_TRACKING)
    #     # self.urukul2_ch1.set_phase_mode(PHASE_MODE_TRACKING)
    #
    #     self.urukul2_ch0.set_phase_mode(PHASE_MODE_ABSOLUTE)
    #     self.urukul2_ch1.set_phase_mode(PHASE_MODE_ABSOLUTE)
    #
    #     self.core.break_realtime()
    #     t = now_mu() + 200_000  # ~2 ms in the future; plenty of slack
    #
    #     self.urukul2_ch0.set(200.0 * MHz, amplitude=(2 * 50 * 10 ** (dBm / 10 - 3)) ** (1 / 2), phase = 0.0, ref_time_mu=t)
    #     self.urukul2_ch1.set(200.0 * MHz, amplitude=(2 * 50 * 10 ** (dBm / 10 - 3)) ** (1 / 2), phase = 0.0, ref_time_mu=t)
    #
    #     # self.urukul2_ch0.set(200.0 * MHz, amplitude=(2 * 50 * 10 ** (dBm / 10 - 3)) ** (1 / 2))
    #     # self.urukul2_ch1.set(200.0 * MHz, amplitude=(2 * 50 * 10 ** (dBm / 10 - 3)) ** (1 / 2))
    #
    #     with parallel:
    #
    #         self.urukul2_ch0.sw.on()
    #         self.urukul2_ch1.sw.on()
    #
    #     # self.urukul2_cpld.io_update.pulse(8 * ns)
    #
    #     delay(100*us)
    #
    #     #
    #     # self.urukul2_ch0.sw.off()
    #     # self.urukul2_ch1.sw.off()
    #
    #     print("code done!")




    def build(self):
       self.setattr_device("core")
       self.setattr_device("urukul2_cpld")
       self.setattr_device("urukul1_cpld")
       self.setattr_device("urukul2_ch0")
       self.setattr_device("urukul1_ch0")

    @kernel
    def run(self):
        self.core.reset()
        self.urukul1_cpld.init()
        self.urukul2_cpld.init()
        self.urukul2_ch0.init()
        self.urukul2_ch0.set_att(float(0))

        self.urukul1_ch0.init()
        self.urukul1_ch0.set_att(float(0))

        delay(1 * ms)

        dBm = -5

        # self.urukul1_ch0.set_phase_mode(PHASE_MODE_TRACKING)
        # self.urukul2_ch0.set_phase_mode(PHASE_MODE_TRACKING)

        self.urukul2_ch0.set_phase_mode(PHASE_MODE_ABSOLUTE)
        self.urukul1_ch0.set_phase_mode(PHASE_MODE_ABSOLUTE)

        self.core.break_realtime()
        t = now_mu() + 200_000  # ~2 ms in the future; plenty of slack

        # self.urukul1_ch0.set(21.0 * MHz, amplitude=(2 * 50 * 10 ** (dBm / 10 - 3)) ** (1 / 2), phase = 0.0, ref_time_mu=t)
        # self.urukul2_ch0.set(21.0 * MHz, amplitude=(2 * 50 * 10 ** (dBm / 10 - 3)) ** (1 / 2), phase = 0.0, ref_time_mu=t)

        self.urukul1_ch0.set(21.0 * MHz, amplitude=(2 * 50 * 10 ** (dBm / 10 - 3)) ** (1 / 2))
        self.urukul2_ch0.set(21.0 * MHz, amplitude=(2 * 50 * 10 ** (dBm / 10 - 3)) ** (1 / 2))

        with parallel:

            self.urukul1_ch0.sw.on()
            self.urukul2_ch0.sw.on()

        delay(500*ns)

        with parallel:

            self.urukul1_ch0.sw.off()
            self.urukul2_ch0.sw.off()

        # self.urukul1_cpld.io_update.pulse(8 * ns)
        # self.urukul2_cpld.io_update.pulse(8 * ns)

        delay(100*us)

        #
        # self.urukul1_ch0.sw.off()
        # self.urukul2_ch0.sw.off()

        print("code done!")



# ### Testing DDSs:
#     def build(self):
#        self.setattr_device("core")
#        self.setattr_device("ttl7")
#
#        self.dds1 = self.get_device("urukul0_ch0")  ### FORT dds
#        self.dds2 = self.get_device("urukul0_ch1")
#
#
#     @kernel
#     def run(self):
#         self.core.reset()
#
#         self.dds1.cpld.init()
#         self.dds1.init()
#         self.dds1.set_att(0.0)
#
#         self.dds1.cpld.set_profile(0)
#
#         delay(10 * us)
#
#         dBm = -5
#         ampl = (2 * 50 * 10 ** (dBm / 10 - 3)) ** (1 / 2)
#         print(ampl)
#         delay(10*ms)
#
#         self.dds1.set(50 * MHz, amplitude=ampl, profile=0)
#         delay(10*us)
#
#         # ### This does not work; no error, but no effect; totally ignored:
#         # self.dds1.set_frequency(5*MHz)
#         # self.dds1.set_amplitude(0.5)
#
#
#         self.ttl7.on()
#         self.dds1.sw.on()
#
#         delay(500 * us)
#         self.dds1.sw.off()
#         self.ttl7.off()
#
#         print("code done!")


 # #### Testing Zotino:
 #    def build(self):
 #        self.setattr_device("core")
 #        self.setattr_device("zotino0")
 #
 #    @kernel
 #    def run(self):
 #        self.core.reset()
 #        self.core.break_realtime()
 #        self.zotino0.init()
 #        delay(10 * ms)
 #
 #        # self.zotino0.write_dac(0, 2.0)
 #        # self.zotino0.write_dac(1, 0.0)
 #        # self.zotino0.write_dac(2, 0.0)
 #        # self.zotino0.write_dac(3, 0.0)
 #        # self.zotino0.write_dac(4, 0.0)
 #        # self.zotino0.write_dac(5, 0.0)
 #        # self.zotino0.load()
 #
 #        for i in range(10000):
 #            self.zotino0.set_dac([2.0], [3])
 #            delay(2*ms)
 #            self.zotino0.set_dac([-2.0], [3])
 #            delay(50*ms)
 #
 #        # self.zotino0.set_dac([2,2,2,2])
 #
 #        delay(10 * ms)
 #
 #        print("code done!")



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






