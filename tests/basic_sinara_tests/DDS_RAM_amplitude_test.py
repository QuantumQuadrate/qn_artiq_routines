# """
# Trying to add a long dwell time between rise and fall RAM modes. Does not work.
#
# Akbar 2025-05-12
#
# """
# from artiq.coredevice.ad9910 import RAM_DEST_ASF, RAM_MODE_RAMPUP
# from artiq.experiment import *
# from artiq.language import us, ns, MHz
# import numpy as np
# import math
#
# class DDS_RAM_amplitude_test(EnvExperiment):
#
#     def build(self):
#         self.setattr_device("core")
#         self.setattr_device("ttl7")
#         self.dds = self.get_device("urukul2_ch0")
#
#     def prepare(self):
#
#         self.ramp_time = 2 * us            # ramp time from first to max amplitude
#         self.N = 100                        # number of amplitude points during the rise and fall time
#         self.amp_low, self.amp_high = 0.0, 1.0       # low and high amplitudes in scale from 0 to 1
#
#         self.step_ticks = int((self.ramp_time / self.N) / (4 * ns))  ### AD9910 can update the RAM every 4 clock cycles, i.e. every 4ns.
#         ### step_ticks=25, for example, means the dds is updated every 25*4ns = 100ns.
#
#         self.dwell_time = 2 * us
#
#         ### Gaussian function. Put this in a python script and plot amp_points (without reversed()) to see the shape of the dds pulse.
#         x_vals = [0.0 + 3.1 * i/(self.N - 1) for i in range(self.N)] ### a list from 0 to 3.1
#         raw    = [math.exp(-0.5 * x*x) for x in x_vals]
#         g_min, g_max = raw[0], raw[-1]
#         norm   = [(r - g_min)/(g_max - g_min) for r in raw]
#
#
#         ### scale into respective ramp according to amplitude.
#         amp_points_fall = [self.amp_low + n * (self.amp_high - self.amp_low) for n in norm]
#         amp_points_rise = list(reversed([self.amp_low + n * (self.amp_high - self.amp_low) for n in norm]))
#
#
#         ### some data conversion needed for RAM
#         amplitudes_rise_arr = np.zeros(len(amp_points_rise), dtype=np.int32)
#         self.dds.amplitude_to_ram(amp_points_rise, amplitudes_rise_arr) ### updates arr according to the amp_points profile
#         self.amplitudes_rise_list = list(amplitudes_rise_arr)
#
#         ### some data conversion needed for RAM
#         amplitudes_fall_arr = np.zeros(len(amp_points_fall), dtype=np.int32)
#         self.dds.amplitude_to_ram(amp_points_fall, amplitudes_fall_arr)  ### updates arr according to the amp_points profile
#         self.amplitudes_fall_list = list(amplitudes_fall_arr)
#
#
#     @kernel
#     def run(self):
#         self.core.reset()
#         self.dds.cpld.init()
#         self.dds.init(blind=True)
#         self.dds.cfg_sw(True)
#         self.core.break_realtime() ### doesn't seem to be necessary
#         self.dds.sw.off()
#         self.dds.set_frequency(10 * MHz)   # the frequency we want
#         self.dds.set_att(0.0)
#
#         self.run_ram(self.step_ticks)
#
#     @kernel
#     def run_ram(self, step_size):
#         # delay(5 * us)
#         self.dds.set_cfr1(ram_enable=0) ### disable RAM mode to write the config
#         self.dds.cpld.io_update.pulse_mu(8) ### pulse the ttl to update and implement settings
#
#
#
#         with parallel:
#             with sequential:
#                 ### Configures the RAM playback engine
#                 self.dds.set_profile_ram(
#                     start=0,
#                     end=self.N - 1,
#                     step=step_size,
#                     profile=0,
#                     mode=RAM_MODE_RAMPUP,
#                 )
#
#                 self.dds.cpld.set_profile(0)
#                 # self.dds.cpld.io_update.pulse_mu(8)
#                 self.dds.write_ram(self.amplitudes_rise_list) ### write the data onto RAM
#
#                 ### Enabling RAM playback, not playing yet.
#                 self.dds.set_cfr1(
#                     internal_profile=0,
#                     ram_enable=1,
#                     ram_destination=RAM_DEST_ASF,
#                 )
#
#                 with parallel:
#                     self.ttl7.on() ### Just for triggering the osc.
#                     self.dds.sw.on()
#                 self.dds.cpld.io_update.pulse_mu(8)  ### This runs the RAM
#
#             with sequential:
#                 delay(self.ramp_time)
#                 ### Configures the RAM playback engine
#                 self.dds.set_profile_ram(
#                     start=0,
#                     end=self.N - 1,
#                     step=step_size,
#                     profile=1,
#                     mode=RAM_MODE_RAMPUP,
#                 )
#
#                 self.dds.cpld.set_profile(1)
#                 self.dds.cpld.io_update.pulse_mu(8)
#                 self.dds.write_ram(self.amplitudes_fall_list)  ### write the data onto RAM
#
#                 ### Enabling RAM playback, not playing yet.
#                 self.dds.set_cfr1(
#                     internal_profile=1,
#                     ram_enable=1,
#                     ram_destination=RAM_DEST_ASF,
#                 )
#                 self.dds.cpld.io_update.pulse_mu(8)  ### This runs the RAM
#
#
#
#
#         # delay(1000 * us)  ### keep the delay as ramp time
#         delay(2*self.ramp_time)  ### keep the delay as ramp time
#
#         ### shutting off
#         self.dds.set_cfr1(ram_enable=0)
#         self.dds.cpld.io_update.pulse_mu(8)
#         self.dds.sw.off()
#         self.ttl7.off()






"""
Based on the code below, but using dds_microwaves to generate smooth MW pulses, send to the mixer, and observed on the scope
from the coupler, after the MW detector. It works and I see a smooth MW pulse on the scope.

Akbar 2025-05-15

"""
from artiq.coredevice.ad9910 import RAM_DEST_ASF, RAM_MODE_RAMPUP
from artiq.experiment import *
from artiq.language import us, ns, MHz
import numpy as np
import math

class DDS_RAM_amplitude_test(EnvExperiment):

    def build(self):
        self.setattr_device("core")
        self.setattr_device("ttl7")
        self.setattr_device("ttl4") ### this is indeed the ttl_microwave_switch
        self.dds = self.get_device("urukul2_ch3")

    def prepare(self):

        self.ramp_time = 20 * us            # ramp time from first to max amplitude
        self.N = 30                        # number of amplitude points during the rise and fall time
        self.amp_low, self.amp_high = 0.0, 0.5       # low and high amplitudes in scale from 0 to 1

        self.step_ticks = int((self.ramp_time / self.N) / (4 * ns))  ### AD9910 can update the RAM every 4 clock cycles, i.e. every 4ns.
        ### step_ticks=25, for example, means the dds is updated every 25*4ns = 100ns.

        self.dwell_time = 0 * us
        self.M = int(self.N /self.ramp_time * self.dwell_time) # Number of dwell points

        ### Gaussian function. Put this in a python script and plot amp_points (without reversed()) to see the shape of the dds pulse.
        x_vals = [0.0 + 3.1 * i/(self.N - 1) for i in range(self.N)] ### a list from 0 to 3.1
        raw    = [math.exp(-0.5 * x*x) for x in x_vals]
        g_min, g_max = raw[0], raw[-1]
        norm   = [(r - g_min)/(g_max - g_min) for r in raw]


        ### scale into respective ramp according to amplitude.
        amp_points_rise = [self.amp_low + n * (self.amp_high - self.amp_low) for n in norm]
        amp_points_fall = list(reversed([self.amp_low + n * (self.amp_high - self.amp_low) for n in norm]))

        ### The full waveform. Start with falling edge!
        amp_points = (
                amp_points_rise +
                [self.amp_high] * self.M +
                amp_points_fall
        )

        print(len(amp_points))


        ### some data conversion needed for RAM
        amplitudes_arr = np.zeros(len(amp_points), dtype=np.int32)
        self.dds.amplitude_to_ram(amp_points, amplitudes_arr) ### updates arr according to the amp_points profile
        self.amplitudes_list = list(amplitudes_arr)

        ### This is calculation of steps based on above parameters
        self.total_points = len(self.amplitudes_list)


    @kernel
    def run(self):
        self.core.reset()
        delay(10 * ms)
        self.ttl4.off()
        delay(10*ms)

        self.dds.cpld.init()
        self.dds.init(blind=True)
        self.dds.init()
        # self.dds.cfg_sw(True) ### this causes a small pulse to be generated about 250us before the desired pulse!
        # self.core.break_realtime() ### doesn't seem to be necessary
        # self.dds.sw.off()

        self.dds.set_att(10.0)
        # self.dds.set(frequency=338.0 * MHz) ### Do not set freq like this. It does not work, though no error!!
        self.dds.set_frequency(338.0 * MHz)  # Use this to set the frequency.
        delay(10 * ms)

        # self.ttl7.pulse(10*us)
        # self.dds.sw.on()
        # delay(10*us)
        # self.dds.sw.off()
        # delay(10*ms)

        self.run_ram(self.step_ticks)

    @kernel
    def run_ram(self, step_size):
        # delay(5 * us)
        self.dds.set_cfr1(ram_enable=0) ### disable RAM mode to write the config
        self.dds.cpld.io_update.pulse_mu(8) ### pulse the ttl to update and implement settings


        ### Configures the RAM playback engine
        self.dds.set_profile_ram(
            start=0,
            end=self.total_points - 1,
            step=step_size,
            profile=0,
            mode=RAM_MODE_RAMPUP,
        )

        self.dds.cpld.set_profile(0)
        # self.dds.cpld.io_update.pulse_mu(8)
        self.dds.write_ram(self.amplitudes_list) ### write the data onto RAM

        ### Enabling RAM playback, not playing yet.
        self.dds.set_cfr1(
            internal_profile=0,
            ram_enable=1,
            ram_destination=RAM_DEST_ASF,
        )

        with parallel:
            self.ttl7.on() ### Just for triggering the osc.
            self.dds.sw.on()
        self.dds.cpld.io_update.pulse_mu(8)  ### This runs the RAM

        delay(50 * us)  ### keep the delay as ramp time
        # delay(self.dwell_time)  ### keep the delay as ramp time

        ### shutting off
        self.dds.set_cfr1(ram_enable=0)
        self.dds.cpld.io_update.pulse_mu(8)
        self.dds.sw.off()
        self.ttl7.off()











# """
# This works well to generate a pulse with smooth rise and fall time. Works well if dwell_time is less than 5us. Other wise, len(amp_points)
# gets too larg (>500) and we get Underflow error.
#
# Akbar 2025-05-12
#
# """
# from artiq.coredevice.ad9910 import RAM_DEST_ASF, RAM_MODE_RAMPUP
# from artiq.experiment import *
# from artiq.language import us, ns, MHz
# import numpy as np
# import math
#
# class DDS_RAM_amplitude_test(EnvExperiment):
#
#     def build(self):
#         self.setattr_device("core")
#         self.setattr_device("ttl7")
#         self.dds = self.get_device("urukul2_ch0")
#
#     def prepare(self):
#
#         self.ramp_time = 2 * us            # ramp time from first to max amplitude
#         self.N = 30                        # number of amplitude points during the rise and fall time
#         self.amp_low, self.amp_high = 0.0, 1.0       # low and high amplitudes in scale from 0 to 1
#
#         self.step_ticks = int((self.ramp_time / self.N) / (4 * ns))  ### AD9910 can update the RAM every 4 clock cycles, i.e. every 4ns.
#         ### step_ticks=25, for example, means the dds is updated every 25*4ns = 100ns.
#
#         self.dwell_time = 24 * us
#         self.M = int(self.N /self.ramp_time * self.dwell_time) # Number of dwell points
#
#         ### Gaussian function. Put this in a python script and plot amp_points (without reversed()) to see the shape of the dds pulse.
#         x_vals = [0.0 + 3.1 * i/(self.N - 1) for i in range(self.N)] ### a list from 0 to 3.1
#         raw    = [math.exp(-0.5 * x*x) for x in x_vals]
#         g_min, g_max = raw[0], raw[-1]
#         norm   = [(r - g_min)/(g_max - g_min) for r in raw]
#
#
#         ### scale into respective ramp according to amplitude.
#         amp_points_rise = [self.amp_low + n * (self.amp_high - self.amp_low) for n in norm]
#         amp_points_fall = list(reversed([self.amp_low + n * (self.amp_high - self.amp_low) for n in norm]))
#
#         ### The full waveform. Start with falling edge!
#         amp_points = (
#                 amp_points_rise +
#                 [self.amp_high] * self.M +
#                 amp_points_fall
#         )
#
#         print(len(amp_points))
#
#
#         ### some data conversion needed for RAM
#         amplitudes_arr = np.zeros(len(amp_points), dtype=np.int32)
#         self.dds.amplitude_to_ram(amp_points, amplitudes_arr) ### updates arr according to the amp_points profile
#         self.amplitudes_list = list(amplitudes_arr)
#
#         ### This is calculation of steps based on above parameters
#         self.total_points = len(self.amplitudes_list)
#
#
#     @kernel
#     def run(self):
#         self.core.reset()
#         self.dds.cpld.init()
#         self.dds.init(blind=True)
#         self.dds.cfg_sw(True) ### this causes a small pulse to be generated about 250us before the desired pulse!
#         # self.core.break_realtime() ### doesn't seem to be necessary
#         self.dds.sw.off()
#         self.dds.set_frequency(100 * MHz)   # the frequency we want
#         self.dds.set_att(20.0)
#
#         self.run_ram(self.step_ticks)
#
#     @kernel
#     def run_ram(self, step_size):
#         # delay(5 * us)
#         self.dds.set_cfr1(ram_enable=0) ### disable RAM mode to write the config
#         self.dds.cpld.io_update.pulse_mu(8) ### pulse the ttl to update and implement settings
#
#
#         ### Configures the RAM playback engine
#         self.dds.set_profile_ram(
#             start=0,
#             end=self.total_points - 1,
#             step=step_size,
#             profile=0,
#             mode=RAM_MODE_RAMPUP,
#         )
#
#         self.dds.cpld.set_profile(0)
#         # self.dds.cpld.io_update.pulse_mu(8)
#         self.dds.write_ram(self.amplitudes_list) ### write the data onto RAM
#         # self.write_ram_safe() ### write the data onto RAM
#
#
#         ### Enabling RAM playback, not playing yet.
#         self.dds.set_cfr1(
#             internal_profile=0,
#             ram_enable=1,
#             ram_destination=RAM_DEST_ASF,
#         )
#
#         with parallel:
#             self.ttl7.on() ### Just for triggering the osc.
#             self.dds.sw.on()
#         self.dds.cpld.io_update.pulse_mu(8)  ### This runs the RAM
#
#         delay(1000 * us)  ### keep the delay as ramp time
#         # delay(self.dwell_time)  ### keep the delay as ramp time
#
#         ### shutting off
#         self.dds.set_cfr1(ram_enable=0)
#         self.dds.cpld.io_update.pulse_mu(8)
#         self.dds.sw.off()
#         self.ttl7.off()









#
# """
# From Akshat. This works well to ramp up a dds amplitude from 0 to 1 within a time that I changed from 2us to 10ms and it worked well.
# Akbar 2025-05-12
#
# """
# from artiq.coredevice.ad9910 import RAM_DEST_ASF, RAM_MODE_RAMPUP
# from artiq.experiment import *
# from artiq.language import us, ns, MHz
# import numpy as np
# import math
#
# class DDS_RAM_amplitude_test(EnvExperiment):
#
#     def build(self):
#         self.setattr_device("core")
#         self.setattr_device("ttl7")
#         self.dds = self.get_device("urukul2_ch0")
#
#     def prepare(self):
#
#         self.ramp_time = 2 * us            # ramp time from first to max amplitude
#         self.N = 100                        # number of amplitude points
#         amp_start, amp_end = 0.0, 1.0       # start and end points in scale from 0 to 1
#
#         ### Gaussian function. Put this in a python script and plot amp_points (without reversed()) to see the shape of the dds pulse.
#         x_vals = [0.0 + 10.0 * i/(self.N - 1) for i in range(self.N)] ### a list from 0 to 10
#         raw    = [math.exp(-0.5 * x*x) for x in x_vals]
#         g_min, g_max = raw[0], raw[-1]
#         norm   = [(r - g_min)/(g_max - g_min) for r in raw]
#
#         ### scale into respective ramp according to amplitude. For some reason, we have to generate the list in reverse.
#         ### The RAM runs the last element of amp_points first, and the first element last.
#         amp_points = list(reversed([amp_start + n*(amp_end - amp_start) for n in norm]))
#
#         ### some data conversion thats needed for RAM
#         amplitudes_arr = np.zeros(self.N, dtype=np.int32)
#         self.dds.amplitude_to_ram(amp_points, amplitudes_arr) ### updates arr according to the amp_points profile
#         self.amplitudes_list = list(amplitudes_arr)
#
#         ### This is calculation of steps based on above parameters
#         self.step_ticks = int((self.ramp_time / self.N) / (4 * ns)) ### AD9910 can update the RAM every 4 clock cycles, i.e. every 4ns.
#         ### step_ticks=25, for example, mean the dds is updated every 25*4ns = 100ns.
#
#     @kernel
#     def run(self):
#         self.core.reset()
#         self.dds.cpld.init()
#         self.dds.init(blind=True)
#         self.dds.cfg_sw(True)
#         # self.core.break_realtime() ### doesn't seem to be necessary
#         self.dds.sw.off()
#         self.dds.set_frequency(10 * MHz)   # the frequency we want
#         self.dds.set_att(0.0)
#
#         self.run_ram(self.step_ticks)
#
#     @kernel
#     def run_ram(self, step_size):
#         delay(5 * us)
#         self.dds.set_cfr1(ram_enable=0) ### disable RAM mode to write the config
#         self.dds.cpld.io_update.pulse_mu(8) ### pulse the ttl to update and implement settings
#
#         ### Configures the RAM playback engine
#         self.dds.set_profile_ram(
#             start=0,
#             end=self.N - 1,
#             step=step_size,
#             profile=0,
#             mode=RAM_MODE_RAMPUP,
#         )
#
#         self.dds.cpld.set_profile(0)
#         self.dds.cpld.io_update.pulse_mu(8)
#         self.dds.write_ram(self.amplitudes_list) ### write the data onto RAM
#
#         ### Enabling RAM playback, not playing yet.
#         self.dds.set_cfr1(
#             internal_profile=0,
#             ram_enable=1,
#             ram_destination=RAM_DEST_ASF,
#         )
#
#         with parallel:
#             self.ttl7.on() ### Just for triggering the osc.
#             self.dds.sw.on()
#         self.dds.cpld.io_update.pulse_mu(8) ### This runs the RAM
#
#         delay(self.ramp_time)  ### keep the delay as ramp time
#
#         ### shutting off
#         self.dds.set_cfr1(ram_enable=0)
#         self.dds.cpld.io_update.pulse_mu(8)
#         self.dds.sw.off()
#         self.ttl7.off()
