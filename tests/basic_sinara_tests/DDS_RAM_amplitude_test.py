"""
Uses dds1 in RAM mode, runs the first half of the RAM to ramp up the signal on dds1, leaves the dds ON for a while, runs the
2nd half of the RAM to ramp down dds1. This works fine. But when I change the settings of dds2
in the middle of the dwell time, that causes dds1 to start the RAM from beginning. This is because both dds are on the same card.
To avoid this issue, we should mask dds1 while we want to change the settings of dds2. Masking is not necessary for just turning
dds2 on/off.

Another point is that if we use profile 7 for RAM, we don't need to specify the profile after the RAM or the profile of other
dds channels. See the 2nd experiment below.

Akbar 2025-06-13
"""

from artiq.experiment import *
from artiq.coredevice.ad9910 import RAM_MODE_RAMPUP, RAM_DEST_ASF
from artiq.coredevice.urukul import CFG_MASK_NU
from artiq.language.types import TInt32
from artiq.language import ms, us, ns, MHz
import numpy as np
import math


class DDS_RAM_amplitude_test(EnvExperiment):

    def build(self):
        self.setattr_device("core")
        self.setattr_device("ttl7")
        self.dds1 = self.get_device("urukul2_ch0")  ### FORT dds
        self.dds2 = self.get_device("urukul2_ch1") ### Cooling DP dds

    def prepare(self):
        self.ramp_time = 0.2 * ms  # ramp time from first to max amplitude
        self.ramp_points = 50  # number of amplitude points during the rise and fall time
        self.amp_low, self.amp_high = 0.0, 0.2  # low and high amplitudes in scale from 0 to 1

        self.step_ticks = int((self.ramp_time / self.ramp_points) / (4 * ns))  ### AD9910 can update the RAM every 4 clock cycles, i.e. every 4ns.
        ### step_ticks=25, for example, means the dds is updated every 25*4ns = 100ns.

        ### Gaussian function. Put this in a python script and plot amp_points (without reversed()) to see the shape of the dds pulse.
        x_vals = [0.0 + 3.1 * i / (self.ramp_points - 1) for i in range(self.ramp_points)]  ### a list from 0 to 3.1
        raw = [math.exp(-0.5 * x * x) for x in x_vals]
        g_min, g_max = raw[0], raw[-1]
        norm = [(r - g_min) / (g_max - g_min) for r in raw]

        ### scale into respective ramp according to amplitude.
        amp_points_rise = [self.amp_low + n * (self.amp_high - self.amp_low) for n in norm]
        amp_points_fall = list(reversed([self.amp_low + n * (self.amp_high - self.amp_low) for n in norm]))

        ### The full waveform
        amp_points = (
                amp_points_rise +
                amp_points_fall
        )

        ### some data conversion needed for RAM
        amplitudes_arr = np.zeros(len(amp_points), dtype=np.int32)
        self.dds1.amplitude_to_ram(amp_points, amplitudes_arr)
        self.ram_data = list(amplitudes_arr)


    @kernel
    def run(self):
        self.core.reset()

        ### dds1 is used in RAM mode
        self.dds1.cpld.init()
        self.dds1.init()
        self.dds1.set_frequency(1.0 * MHz)  # Use this to set the frequency.
        ### Do not set freq like this: self.dds.set(frequency=338.0 * MHz).  It does not work for RAM, though no error!!
        self.dds1.set_att(0.0)
        self.dds1.sw.off()

        ### dds2 is used in non-RAM mode
        self.dds2.cpld.init()
        self.dds2.init()
        self.dds2.set(frequency=1.0 * MHz, amplitude=0.1)
        self.dds2.set_att(0.0)
        self.dds2.sw.off()

        delay(10 * us)

        self.dds1.set_cfr1(ram_enable=0)  ### disable RAM mode to write the config
        self.dds1.cpld.io_update.pulse_mu(8)  ### pulse the ttl to update and implement settings

        ### Configures the RAM playback engine
        self.dds1.set_profile_ram(
            start=0,
            end=len(self.ram_data) - 1,
            step=self.step_ticks,
            profile=7,
            mode=RAM_MODE_RAMPUP,
        )

        ### write the data onto RAM
        # self.dds1.cpld.set_profile(1)
        self.dds1.cpld.io_update.pulse_mu(8)
        self.dds1.write_ram(self.ram_data)

        delay(0.1 * ms)

        self.ttl7.on()


        ### Configure the RAM to playback the first half
        self.dds1.set_profile_ram(
            start=0,
            end=len(self.ram_data)//2 -1,
            step=self.step_ticks,
            profile=7,
            mode=RAM_MODE_RAMPUP)
        # self.dds1.cpld.set_profile(1)
        self.dds1.cpld.io_update.pulse_mu(8)

        self.core.break_realtime()

        ### Enabling RAM playback, not playing yet.
        self.dds1.set_cfr1(ram_enable=1,
                          ram_destination=RAM_DEST_ASF)

        ### Running the RAM
        self.dds1.sw.on()
        self.dds1.cpld.io_update.pulse_mu(8)
        delay(0.5*ms)  # Leave at-least enough time to cover up the first ram time, plus any extra time we want.

        ### We don't need to mask dds1 to turn on or off dds2. But we need to mask it when changing the settings of dds2.
        ### Otherwise, if we disable the two mask lines below, dds1 RAM restarts when we change dds2 setting.
        self.dds2.sw.on()
        delay(0.2*ms)
        self.dds2.sw.off()
        delay(0.2 * ms)
        self.dds2.sw.on()

        delay(0.5 * ms)
        self.dds1.cpld.cfg_write(self.dds1.cpld.cfg_reg | 1 << CFG_MASK_NU + 0)  # Mask the DDS channel which is on RAM mode
        self.dds2.set(frequency=1.0 * MHz, amplitude=0.1)
        self.dds1.cpld.cfg_write(self.dds1.cpld.cfg_reg & ~(1 << CFG_MASK_NU + 0))  # Unmask the DDS channel which in on RAM mode
        delay(0.5*ms)


        ### Configure the RAM to playback the second half
        self.dds1.set_profile_ram(
            start=len(self.ram_data)//2,
            end=len(self.ram_data)-1,
            step=self.step_ticks,
            profile=7,
            mode=RAM_MODE_RAMPUP)
        # self.dds1.cpld.set_profile(1)
        self.dds1.cpld.io_update.pulse_mu(8)

        self.dds2.sw.off()

        self.ttl7.off()

        delay(0.5 * ms)

        self.cfr_reset(self.dds1) # Exit RAM Mode

        ### to show that we can change the dds settings after RAM. Note that we have to set the frequency and amplitude
        ### after exiting the RAM (once is enough). Otherwise, the dds will not turn on.
        self.dds1.set(frequency=1.0 * MHz, amplitude=0.1)
        self.dds1.sw.on()
        delay(100*us)
        self.dds1.sw.off()
        delay(100*us)
        self.dds1.sw.on()
        delay(100 * us)
        self.dds1.sw.off()

        delay(0.4 * ms)

        self.dds1.sw.off()

        # self.dds1.sw.off()

    @kernel
    def cfr_reset(self, dds):
        dds.set_cfr1(ram_enable=0)
        dds.cpld.io_update.pulse_mu(8)





# """
# Uses dds1 in RAM mode, runs the first half of the RAM to ramp up the signal on dds1, leaves the dds ON for a while, runs the
# 2nd half of the RAM to ramp down dds1. This works fine as long as dds2 is not used. but when I change the settings of dds2
# in the middle of the dwell time, that causes dds1 to start the RAM from beginning!!
#
# You can use this example for how to run RAM in a profile other than 7.
#
# Akbar 2025-06-06
# """
#
# from artiq.experiment import *
# from artiq.coredevice.ad9910 import RAM_MODE_RAMPUP, RAM_DEST_ASF
# from artiq.language.types import TInt32
# from artiq.language import ms, us, ns, MHz
# import numpy as np
# import math
#
#
# class DDS_RAM_amplitude_test(EnvExperiment):
#
#     def build(self):
#         self.setattr_device("core")
#         self.setattr_device("ttl7")
#         self.dds1 = self.get_device("urukul2_ch0")  ### FORT dds
#         self.dds2 = self.get_device("urukul2_ch1") ### Cooling DP dds
#
#     def prepare(self):
#         self.ramp_time = 5 * ms  # ramp time from first to max amplitude
#         self.ramp_points = 50  # number of amplitude points during the rise and fall time
#         self.amp_low, self.amp_high = 0.0, 0.2  # low and high amplitudes in scale from 0 to 1
#
#         self.step_ticks = int((self.ramp_time / self.ramp_points) / (4 * ns))  ### AD9910 can update the RAM every 4 clock cycles, i.e. every 4ns.
#         ### step_ticks=25, for example, means the dds is updated every 25*4ns = 100ns.
#
#         ### Gaussian function. Put this in a python script and plot amp_points (without reversed()) to see the shape of the dds pulse.
#         x_vals = [0.0 + 3.1 * i / (self.ramp_points - 1) for i in range(self.ramp_points)]  ### a list from 0 to 3.1
#         raw = [math.exp(-0.5 * x * x) for x in x_vals]
#         g_min, g_max = raw[0], raw[-1]
#         norm = [(r - g_min) / (g_max - g_min) for r in raw]
#
#         ### scale into respective ramp according to amplitude.
#         amp_points_rise = [self.amp_low + n * (self.amp_high - self.amp_low) for n in norm]
#         amp_points_fall = list(reversed([self.amp_low + n * (self.amp_high - self.amp_low) for n in norm]))
#
#         ### The full waveform
#         amp_points = (
#                 amp_points_rise +
#                 amp_points_fall
#         )
#
#         ### some data conversion needed for RAM
#         amplitudes_arr = np.zeros(len(amp_points), dtype=np.int32)
#         self.dds1.amplitude_to_ram(amp_points, amplitudes_arr)
#         self.ram_data = list(amplitudes_arr)
#
#
#     @kernel
#     def run(self):
#         self.core.reset()
#
#         ### dds1 is used in RAM mode
#         self.dds1.cpld.init()
#         self.dds1.init()
#         self.dds1.set_frequency(1.0 * MHz)  # Use this to set the frequency.
#         ### Do not set freq like this: self.dds.set(frequency=338.0 * MHz).  It does not work for RAM, though no error!!
#         self.dds1.set_att(0.0)
#         self.dds1.sw.off()
#
#         ### dds2 is used in non-RAM mode
#         self.dds2.cpld.init()
#         self.dds2.init()
#         self.dds2.set(frequency=1 * MHz, amplitude=0.1)
#         self.dds2.set_att(0.0)
#         self.dds2.sw.off()
#
#         delay(10 * us)
#
#         self.dds1.set_cfr1(ram_enable=0)  ### disable RAM mode to write the config
#         self.dds1.cpld.io_update.pulse_mu(8)  ### pulse the ttl to update and implement settings
#
#         ### Configures the RAM playback engine
#         self.dds1.set_profile_ram(
#             start=0,
#             end=len(self.ram_data) - 1,
#             step=self.step_ticks,
#             profile=1,
#             mode=RAM_MODE_RAMPUP,
#         )
#
#         ### write the data onto RAM
#         self.dds1.cpld.set_profile(1)
#         self.dds1.cpld.io_update.pulse_mu(8)
#         self.dds1.write_ram(self.ram_data)
#
#         delay(1 * ms)
#
#         self.ttl7.on()
#
#
#         ### Configure the RAM to playback the first half
#         self.dds1.set_profile_ram(
#             start=0,
#             end=len(self.ram_data)//2 -1,
#             step=self.step_ticks,
#             profile=1,
#             mode=RAM_MODE_RAMPUP)
#         self.dds1.cpld.set_profile(1)
#         self.dds1.cpld.io_update.pulse_mu(8)
#
#         self.core.break_realtime()
#
#         ### Enabling RAM playback, not playing yet.
#         self.dds1.set_cfr1(ram_enable=1,
#                           ram_destination=RAM_DEST_ASF)
#
#         ### Running the RAM
#         self.dds1.sw.on()
#         self.dds1.cpld.io_update.pulse_mu(8)
#
#         self.dds2.sw.on()
#
#         delay(5 * ms)   # Leave at-least enough time to cover up the first ram time, plus any extra time we want.
#         self.dds2.set(frequency=1 * MHz, amplitude=0.5, profile=0) ### this messes up dds1 and runs the RAM again!!
#         delay(5*ms)
#
#
#         ### Configure the RAM to playback the second half
#         self.dds1.set_profile_ram(
#             start=len(self.ram_data)//2,
#             end=len(self.ram_data)-1,
#             step=self.step_ticks,
#             profile=1,
#             mode=RAM_MODE_RAMPUP)
#         self.dds1.cpld.set_profile(1)
#         self.dds1.cpld.io_update.pulse_mu(8)
#
#         self.dds2.sw.off()
#
#         self.ttl7.off()
#
#         delay(5 * ms)
#
#         self.cfr_reset(self.dds1) # Exit RAM Mode
#
#         # ### to show that we can change the dds settings after RAM
#         # self.dds1.set(frequency=245 * MHz, amplitude=0.1, profile=1)
#         # self.dds1.sw.on()
#
#         delay(4 * ms)
#
#         self.dds1.sw.off()
#
#         # self.dds1.sw.off()
#
#     @kernel
#     def cfr_reset(self, dds):
#         dds.set_cfr1(ram_enable=0)
#         dds.cpld.io_update.pulse_mu(8)