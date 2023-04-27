# this syntax denotes the ARTIQ-type equivalent of the return and must be used in order to return something
def volts_to_power_test(V) -> TFloat:
    P = V*0.67 + 0.5
    return P

# can import in external exp - works
@kernel
def measure(sampler):
    smp = [0.0]*8
    sampler.sample(smp)

# can import in external exp - works
@kernel
def set_dds(dds_channel, freq, amp):
    dds_channel.set(freq, amplitude=amp)

@kernel
def dummy_conversion(x: TFloat) -> TFloat:
    return x

# class ServoChannels:
#
#     # @rpc(flags={'async'})
#     def __init__(self,dds_channels=[],
#                  sampler=None,
#                  sampler_channels=[],
#                  conversions=[],
#                  p_coeffs=[],
#                  setpoints=[]):
#         """
#         A class for adjusting DDS channel powers in response to a Sampler measurement.
#
#         Exactly one instance of this class should be used per sampler.
#
#         :param dds_channels: a list of Urukul channels with len(dds_channels) <= 8
#         :param sampler: a Sampler device
#         :param sampler_channels: a list of the channel indices 0-7 for the sampler, where the
#             channel index at each list position corresponds to the same servo loop as the
#             dds channel at the same position in the dds_channels list.
#         :param conversions: a list of callables which return floats. one per Urukul channel.
#             these convert a measured voltage into a more useful unit, e.g. a scaled voltage
#             or optical power in mW for a given detector. the output unit should match the
#             units of the setpoints.
#         :param setpoints: a list of float values
#         """
#         self.dds_channels = dds_channels
#         # dds_settings = [ch.get() for ch in self.dds_channels]
#         # self.dds_freqs = [settings[0] for settings in dds_settings]
#         # self.dds_amplitudes = [settings[2] for settings in dds_settings]
#         self.sampler = sampler
#         self.sampler_channels = sampler_channels
#         self.n_channels = len(self.dds_channels)
#         self.conversions = conversions # convert what we measure, i.e. V to e.g. mW
#         self.setpoints = setpoints # the values we want in the sensible units, e.g. mW
#         self.buffer = [0.0]*8 # measurement buffer to pass to Sampler
#         self.signals = [6.34]*self.n_channels # will store the measurements
#         self.p_coeffs = p_coeffs
#         self.errors = [0.0]*self.n_channels
#
#     def setup(self):
#         pass
#
#     @kernel
#     def measure(self):
#         self.sampler.sample(self.buffer)
#         # self.signals = [self.buffer[i] for i in self.sampler_channels]
#         self.signals = [(7/0.9)*ampl for ampl in self.dds_amplitudes]
#
#
#     def update_dds_amplitudes(self):
#         """
#         this updates the class variable that keeps a record -
#         it does not set the output of the DDS channels
#         """
#         self.dds_amplitudes = [ch.get()[2] for ch in self.dds_channels]
#
#
#     @rpc(flags={"async"})
#     def print_error(self):
#         print(self.errors)
#
#
#     @kernel
#     def feedback(self):
#
#         pass
#
#         # for i in range(10): # to-do specify tolerance we want to reach?
#         #
#         #     self.measure()
#         #     self.errors = [sp-conv(sig) for sp,conv,sig in zip(self.setpoints,self.conversions,self.signals)]
#         #     self.print_error()
#         #
#         #     for ch, freq, err, p in zip(self.dds_channels, self.dds_freqs, self.errors, self.p_coeffs):
#         #         ampl = ch.get()[2] + p*err
#         #         ampl = ampl if ampl < 1 else 0.9 # don't try to drive the DDS beyond its limit
#         #         ch.set(freq, amplitude=ampl)
#         #     self.update_dds_amplitudes()