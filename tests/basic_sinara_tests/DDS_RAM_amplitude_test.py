
# ############ from Akshat. Generates a pulse, but does not ramp smoothly.
#
# from artiq.experiment import EnvExperiment, kernel, delay
# from artiq.language.units import us
# from artiq.coredevice.ad9910 import RAM_MODE_CONT_RAMPUP, RAM_DEST_ASF
# from numpy import int32
#
# class DDS_RAM_amplitude_test(EnvExperiment):
#     def build(self):
#         self.setattr_device("core")
#         self.setattr_device("urukul2_cpld")
#         self.setattr_device("urukul2_ch0")
#
#     def prepare(self):
#         self.N         = 150          # number of amplitude steps
#         self.freq      = 10e6        # frequency set
#         self.ramp_time = 20 * us      # ramp time
#
#     @kernel
#     def run(self):
#         core = self.core
#         cpld = self.urukul2_cpld
#         dds  = self.urukul2_ch0
#
#         core.reset()
#         cpld.init()
#         dds.init(blind=True)
#         dds.cfg_sw(True)
#         start_amp = 0.10
#         end_amp   = 1.00
#         amp_list = [
#             start_amp + (end_amp - start_amp) * i / (self.N - 1)
#             for i in range(self.N)
#         ]
#
#         ram_data = [int32(0)] * self.N
#
#         dds.set_cfr1(ram_enable=0)    # disable OSK
#         cpld.io_update.pulse_mu(8)
#         dds.set_profile_ram(
#             start=0,
#             end=self.N - 1,
#             step=1,
#             profile=0,
#             mode=RAM_MODE_CONT_RAMPUP
#         )
#         cpld.io_update.pulse_mu(8)
#
#         dds.amplitude_to_ram(amp_list, ram_data)
#         dds.write_ram(ram_data)
#
#         dds.set(frequency=self.freq)
#         dds.set_cfr1(
#             ram_enable=1,
#             ram_destination=RAM_DEST_ASF
#         )
#         cpld.io_update.pulse_mu(8)
#
#         delay(self.ramp_time)
#
#         dds.set_cfr1(ram_enable=0)
#         cpld.io_update.pulse_mu(8)
#         dds.cfg_sw(False)
#
#         delay(self.ramp_time)
#         self.core.reset()
#         self.urukul2_cpld.init()
#         self.urukul2_ch0.init()
#         self.urukul2_ch0.set_att(0.0)






# ########## from artiq manual. Does not work.
# from artiq.experiment import *
# from artiq.coredevice.ad9910 import RAM_MODE_CONT_RAMPUP
#
# class DDS_RAM_amplitude_test(EnvExperiment):
#     def prepare(self):
#         self.amp = [0.0, 0.0, 0.0, 0.7, 0.0, 0.7, 0.7]  # Reversed Order
#         self.asf_ram = [0] * len(self.amp)
#
#     def build(self):
#         self.setattr_device("core")
#         self.setattr_device("urukul2_cpld")
#         self.setattr_device("urukul2_ch0")
#
#     @kernel
#     def init_dds(self, dds):
#         self.core.break_realtime()
#         dds.init()
#         dds.set_att(6. * dB)
#         dds.cfg_sw(True)
#
#     @kernel
#     def configure_ram_mode(self, dds):
#         self.core.break_realtime()
#         dds.set_cfr1(ram_enable=0)
#         self.cpld.io_update.pulse_mu(8)
#         self.cpld.set_profile(0)  # Enable the corresponding RAM profile
#         # Profile 0 is the default
#         dds.set_profile_ram(start=0, end=len(self.asf_ram) - 1,
#                             step=250, profile=0, mode=RAM_MODE_CONT_RAMPUP)
#         self.cpld.io_update.pulse_mu(8)
#         dds.amplitude_to_ram(self.amp, self.asf_ram)
#         dds.write_ram(self.asf_ram)
#         self.core.break_realtime()
#         dds.set(frequency=5 * MHz, ram_destination=RAM_DEST_ASF)
#         # Pass osk_enable=1 to set_cfr1() if it is not an amplitude RAM
#         dds.set_cfr1(ram_enable=1, ram_destination=RAM_DEST_ASF)
#         self.cpld.io_update.pulse_mu(8)
#
#
#     @kernel
#     def run(self):
#         self.core.reset()
#         self.core.break_realtime()
#         self.urukul2_cpld.init()
#         self.urukul2_ch0.init()
#         self.configure_ram_mode(self.urukul2_ch0)



















########### from gpt. Gives no error, but no output on scope either:
# from artiq.experiment import *
# from artiq.coredevice.ad9910 import AD9910, RAM_MODE_CONT_RAMPUP, RAM_DEST_ASF
#
# class DDS_RAM_amplitude_test(EnvExperiment):
#     def build(self):
#         self.setattr_device("core")
#         self.setattr_device("urukul2_cpld")
#         self.setattr_device("urukul2_ch0")
#
#     @kernel
#     def run(self):
#         self.core.reset()
#         self.urukul2_cpld.init()
#         self.urukul2_ch0.init()
#         self.urukul2_ch0.set_att(0.0)
#
#         # Set the RAM profile with Gaussian pulse (example data)
#         ram_data = [int(1000 * (2.718281828459045 ** (-((i - 50) ** 2) / 100))) for i in range(100)]
#
#         # Set up RAM profile (using the mode constant from the manual example)
#         self.urukul2_ch0.set_profile_ram(
#             start=0,
#             end=len(ram_data) - 1,
#             step=1,
#             profile=0,
#             mode=RAM_MODE_CONT_RAMPUP  # Use the constant from the manual
#         )
#
#         # Load data into RAM
#         self.urukul2_ch0.write_ram(ram_data)  # Use write_ram to load the data
#
#         # Set frequency and RAM destination (from example)
#         self.urukul2_ch0.set(frequency=50 * MHz, ram_destination=RAM_DEST_ASF)
#
#         # Enable RAM (if it's not an amplitude RAM, use the osk_enable flag)
#         self.urukul2_ch0.set_cfr1(ram_enable=1, ram_destination=RAM_DEST_ASF)
#
#         # Update the CPLD
#         self.urukul2_cpld.io_update.pulse_mu(8)
#
#         # Play the RAM profile
#         self.urukul2_ch0.sw.on()
#         delay(5 * ms)
#         self.urukul2_ch0.sw.off()
#
#         delay(10*ms)
#         self.core.reset()
#         self.urukul2_cpld.init()
#         self.urukul2_ch0.init()
#         self.urukul2_ch0.set_att(0.0)
