#this experiment shows that we used RAM in the starting and then successfully closed the ram and used DDS outside RAM
#to use the DDS on different frequency, dds, or phase
from artiq.coredevice.ad9910 import RAM_DEST_ASF, RAM_MODE_RAMPUP, PHASE_MODE_ABSOLUTE, PHASE_MODE_CONTINUOUS, PHASE_MODE_TRACKING
from artiq.experiment import EnvExperiment, kernel
from artiq.language import us, ns, MHz, ms
import numpy as np
import math
class DDS_RAM_amplitude_test(EnvExperiment):

    def build(self):
        self.setattr_device("core")
        self.setattr_device("ttl7")
        self.dds = self.get_device("urukul0_ch0")
        self.dds2 = self.get_device("urukul0_ch1")
        self.setattr_device("urukul0_cpld")

    def prepare(self):

        ramp_time = 20 * us                  # ramp time from first to last amplitude
        N = 100                              # number of amplitude points
        amp_start, amp_end = 1.0, 0.0       # The ending amplitude write first, starting amplitude second one

        # gaussian function
        x_vals = [-3.0 + 3.0 * i/(N - 1) for i in range(N)]
        raw    = [math.exp(-0.5 * x*x)       for x in x_vals]
        g_min, g_max = raw[0], raw[-1]
        norm   = [(r - g_min)/(g_max - g_min) for r in raw]

        # scale into respective ramp according to amplitude
        amp_points = [amp_start + n*(amp_end - amp_start) for n in norm]

        # some data conversion thats needed for RAM
        arr = np.zeros(N, dtype=np.int32)
        self.dds.amplitude_to_ram(amp_points, arr)
        self.data1 = list(arr)

        # This is  calculation of steps based on above parameters
        step_ticks = int((ramp_time / N) / (4 * ns))
        self.steps = [step_ticks]

    @kernel
    def run(self):
        self.core.reset()
        self.dds.cpld.init()
        self.dds.init()
        self.core.break_realtime()
        self.dds.set_frequency(245 * MHz)  # the frequency we want
        self.dds.set_att(0.0)
        self.dds2.set_att(0.0)

        for step in self.steps:
            self.run_ram(step)

        self.ram_mode_exit()

        self.dds.set(frequency=245 * MHz, amplitude=0.5,
                     profile=0,
                     phase=0.5)
        # self.dds.cpld.io_update.pulse_mu(8)
        self.dds.sw.on()

        delay(20 * us)


        self.dds.sw.off()
        self.ttl7.off()
        self.dds2.sw.off()

    @kernel
    def run_ram(self, timestep_mu):
        # delay(5 * us)
        self.dds.set_cfr1(ram_enable=0)
        self.dds.cpld.io_update.pulse_mu(8)
        self.dds.set_profile_ram(
            start=0,
            end=len(self.data1) - 1,
            step=timestep_mu,
            profile=0,
            mode=RAM_MODE_RAMPUP,
        )
        self.dds.cpld.set_profile(0)
        self.dds.cpld.io_update.pulse_mu(8)
        self.dds.write_ram(self.data1)
        self.dds.set_cfr1(
            # internal_profile=0,
            ram_enable=1,
            ram_destination=RAM_DEST_ASF,
        )
        with parallel:
            self.ttl7.on()
            self.dds.sw.on()
        self.dds.cpld.io_update.pulse_mu(8)


        delay(10 * us)  # keep the time as ramp time

        # shutting off


    @kernel
    def cfr_reset(self, dds):
        delay(10 * us)
        self.dds.set_cfr1(
            # internal_profile=0,
            ram_enable=0)
        self.dds.cpld.io_update.pulse_mu(8)  # doesn't matter which dds we set this to

    @kernel
    def ram_mode_exit(self):
        self.cfr_reset(self.dds)  # Exit RAM Mode






#
# """
# This runs the RAM in two separate instances with a delay in between so we can split the RAM, ram up first,
# have an arbitrary delay, then ramp down.
# Initially I got RecursionError: maximum recursion depth exceeded, but the same code worked later!
#
# Akbar 2025-06-04
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
#         self.dds1 = self.get_device("urukul0_ch0")  ### FORT dds
#         self.dds2 = self.get_device("urukul2_ch0")
#
#     def prepare(self):
#         self.ramp_time = 2 * ms  # ramp time from first to max amplitude
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
#         ### The full waveform. Start with falling edge!
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
#         self.dds1.set_frequency(245.0 * MHz)  # Use this to set the frequency.
#         ### Do not set freq like this: self.dds.set(frequency=338.0 * MHz).  It does not work for RAM, though no error!!
#         self.dds1.set_att(0.0)
#         self.dds1.sw.off()
#
#         ### dds2 is used in non-RAM mode
#         self.dds2.cpld.init()
#         self.dds2.init()
#         self.dds2.set(frequency=10 * MHz, amplitude=0.1, profile=0)
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
#             profile=0,
#             mode=RAM_MODE_RAMPUP,
#         )
#
#         ### write the data onto RAM
#         self.dds1.cpld.set_profile(0)
#         self.dds1.cpld.io_update.pulse_mu(8)
#         self.dds1.write_ram(self.ram_data)
#
#         delay(1 * ms)
#
#         self.ttl7.on()
#
#
#
#         ### Configure the RAM to playback the first half
#         self.dds1.set_profile_ram(
#             start=0,
#             end=len(self.ram_data)//2 -1,
#             step=self.step_ticks,
#             profile=0,
#             mode=RAM_MODE_RAMPUP)
#         self.dds1.cpld.set_profile(0)
#         self.dds1.cpld.io_update.pulse_mu(8)
#
#         self.core.break_realtime()
#
#         ### Enabling RAM playback, not playing yet.
#         self.dds1.set_cfr1(internal_profile=0,
#                           ram_enable=1,
#                           ram_destination=RAM_DEST_ASF)
#
#         ### Running the RAM
#         self.dds1.sw.on()
#         self.dds1.cpld.io_update.pulse_mu(8)
#
#         self.dds2.sw.on()
#
#
#
#         delay(8 * ms)   # Leave at-least enough time to cover up the first ram time, plus any extra time we want.
#
#
#
#         ### Configure the RAM to playback the second half
#         self.dds1.set_profile_ram(
#             start=len(self.ram_data)//2,
#             end=len(self.ram_data)-1,
#             step=self.step_ticks,
#             profile=0,
#             mode=RAM_MODE_RAMPUP)
#         self.dds1.cpld.set_profile(0)
#         self.dds1.cpld.io_update.pulse_mu(8)
#
#         self.dds2.sw.off()
#
#
#
#         self.ttl7.off()
#
#         delay(10 * us)
#
#         # self.cfr_reset(self.dds1) # Exit RAM Mode
#
#         # self.dds1.sw.off()
#
#     @kernel
#     def cfr_reset(self, dds):
#         dds.set_cfr1(internal_profile=0, ram_enable=0)
#         dds.cpld.io_update.pulse_mu(8)











#
# """
# This experiment runs gaussain amplitude sweep with RAM profile on one channel (dds1) while using another channel on the same
# card (dds2) in Non-RAM mode. This generates smooth MW pulses, sent to the mixer, and we can observe on the scope
# from the coupler, after the MW detector. It works and I see a smooth MW pulse on the scope.
#
# Akbar 2025-06-04
# """
#
#
# from artiq.coredevice.ad9910 import RAM_DEST_ASF, RAM_MODE_RAMPUP
# from artiq.experiment import EnvExperiment, kernel, delay, parallel, delay_mu
# from artiq.language import ms, us, ns, MHz
# import numpy as np
# import math
#
# class DDS_RAM_amplitude_test(EnvExperiment):
#
#     def build(self):
#         self.setattr_device("core")
#         self.setattr_device("ttl7")
#         self.setattr_device("ttl4") ### ttl_microwave_switch
#         self.dds1 = self.get_device("urukul2_ch3") ### microwave dds
#         self.dds2 = self.get_device("urukul2_ch0")
#
#     def prepare(self):
#         self.ramp_time = 5 * us  # ramp time from first to max amplitude
#         self.ramp_points = 30  # number of amplitude points during the rise and fall time
#         self.amp_low, self.amp_high = 0.0, 0.16  # low and high amplitudes in scale from 0 to 1
#
#         self.step_ticks = int((self.ramp_time / self.ramp_points) / (4 * ns))  ### AD9910 can update the RAM every 4 clock cycles, i.e. every 4ns.
#         ### step_ticks=25, for example, means the dds is updated every 25*4ns = 100ns.
#
#         self.dwell_time = 5 * us
#         self.dwell_points = int(self.ramp_points / self.ramp_time * self.dwell_time)  # Number of dwell points
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
#         ### The full waveform. Start with falling edge!
#         amp_points = (
#                 amp_points_rise +
#                 [self.amp_high] * self.dwell_points +
#                 amp_points_fall
#         )
#
#         ### some data conversion needed for RAM
#         amplitudes_arr = np.zeros(len(amp_points), dtype=np.int32)
#         self.dds1.amplitude_to_ram(amp_points, amplitudes_arr)
#         self.ram_data = list(amplitudes_arr)
#
#     @kernel
#     def run(self):
#         self.core.reset()
#         delay(1*ms)
#
#         ### dds1 is used in RAM mode
#         self.dds1.cpld.init()
#         self.dds1.init() ### This turns on the dds channel for about 60 ms.
#         self.dds1.set_frequency(334.0 * MHz) # Use this to set the frequency.
#         ### Do not set freq like this: self.dds.set(frequency=338.0 * MHz).  It does not work for RAM, though no error!!
#         self.dds1.set_att(0.0)
#         self.dds1.sw.off()
#
#         ### dds2 is used in non-RAM mode
#         self.dds2.cpld.init()
#         self.dds2.init()
#         self.dds2.set(frequency=1 * MHz, amplitude=0.1, profile=0)
#         self.dds2.set_att(0.0)
#         self.dds2.sw.off()
#
#         self.ttl4.off() ### turning on the microwave switch to let MW pass
#         delay(10*us)
#
#         self.play_single_ram()
#
#         self.cfr_reset(self.dds1) # Exit RAM Mode
#
#         self.ttl4.on()
#
#     @kernel
#     def play_single_ram(self):
#
#         self.dds1.set_cfr1(ram_enable=0) ### disable RAM mode to write the config
#         self.dds1.cpld.io_update.pulse_mu(8) ### pulse the ttl to update and implement settings
#
#         ### Configures the RAM playback engine
#         self.dds1.set_profile_ram(
#             start=0,
#             end=len(self.ram_data) - 1,
#             step=self.step_ticks,
#             profile=0,
#             mode=RAM_MODE_RAMPUP,
#         )
#
#         ### write the data onto RAM
#         self.dds1.cpld.set_profile(0)
#         self.dds1.cpld.io_update.pulse_mu(8)
#         self.dds1.write_ram(self.ram_data)
#
#         ### Enabling RAM playback, not playing yet.
#         self.dds1.set_cfr1(internal_profile=0,
#                           ram_enable=1,
#                           ram_destination=RAM_DEST_ASF)
#
#         for _ in range(3):
#             self.dds1.sw.on()
#             self.dds1.cpld.io_update.pulse_mu(8) ### This runs the RAM
#
#             delay(3 * us)
#
#             self.dds2.sw.on()
#             delay(20 * us)
#             self.dds2.sw.off()
#             delay(20 * us)
#             ### the total delay in the loop should be larger than the RAM pulse length. Otherwise, the pulses get cut.
#
#         # for _ in range(3):
#         #     self.dds1.sw.on()
#         #     self.dds1.cpld.io_update.pulse_mu(8) ### This runs the RAM
#         #     delay(2*self.ramp_time + self.dwell_time) ### the delay should be the pulse length. Otherwise, the next
#         #     ### loop starts before the first pulse is done. If delay is not enough, the pulse shapes will be affected.
#
#     @kernel
#     def cfr_reset(self, dds):
#         dds.set_cfr1(internal_profile=0, ram_enable=0)
#         dds.cpld.io_update.pulse_mu(8)