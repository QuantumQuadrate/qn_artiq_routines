"""
Read the frequency output by the Rigol DS1022Z using SA44B Signal Hound function generator

This is useful if we are controlling the Rigol frequency with a voltage input via the modulation port,
as the actual output is not linear in the modulation voltage. This lets us measure the actual output.
"""

from artiq.experiment import *
import logging
import numpy as np
from scipy.signal import get_window

import pyvisa as visa
import os, sys

cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"repository\\qn_artiq_routines")
from third_party.signal_hound.sadevice.sa_api import *


class ReadRigolFrequencies(EnvExperiment):

    def build(self):
        self.setattr_device("core")
        self.setattr_device("zotino0")

        self.setattr_argument("modulation_frequency_sequence", StringValue(
            'np.array([])'))
        self.setattr_argument("Rigol_carrier_frequency", NumberValue(65*kHz, unit='kHz', ndecimals=1))
        self.setattr_argument("Rigol_FM_deviation", NumberValue(65*kHz, unit='kHz', ndecimals=1))

    def prepare(self):
        self.scan_sequence1 = eval(self.modulation_frequency_sequence)

        # make sure we aren't going to generate a voltage outside of the -5 to 5V range of the Rigol input
        assert min(self.scan_sequence1) >= self.Rigol_carrier_frequency - self.Rigol_FM_deviation, \
            "f_modulation lower bound should be > carrier_frequency - FM_deviation!"
        assert max(self.scan_sequence1) <= self.Rigol_carrier_frequency + self.Rigol_FM_deviation, \
            "f_modulation upper bound should be < carrier_frequency + FM_deviation!"

        # this list is only used for the asserts below
        V_modulation_list = (self.scan_sequence1 - self.Rigol_carrier_frequency) * 5 / self.Rigol_FM_deviation

        # this should never be triggered because the asserts above should catch this,
        # but it serves as redundancy in case someone changes the frequency modulation definitions
        assert min(V_modulation_list) >= -5, f"{min(V_modulation_list)} is an invalid voltage for Rigol"
        assert max(V_modulation_list) <= 5, f"{max(V_modulation_list)} is an invalid voltage for Rigol"

        self.actual_frequency_list = np.zeros(len(self.scan_sequence1))

        # the instrument name as found in NI MAX or in the list above
        NAME = 'USB0::0x1AB1::0x0642::DG1ZA252402452::INSTR'

        # # create a VISA reference just to make sure it is closed before starting the experiment
        # rm = visa.ResourceManager()
        # funcgen = rm.open_resource(NAME, timeout=20, chunk_size=1024000)
        # funcgen.close()
        # rm.close()

    # @rpc(flags={"async"})
    def query_Rigol_frequency(self) -> TFloat:
        """
        Get the frequency from a connected Rigol function generator

        :return TFloat: the frequency of the Rigol DG1022Z in Hz
        """

        # querying the device with the Visa connection doesn't work because
        # the modulation input seems to be ignored if we have the USB plugged in
        # # the instrument name as found in NI MAX or in the list above
        # NAME = 'USB0::0x1AB1::0x0642::DG1ZA252402452::INSTR'
        #
        # # create a VISA reference to the scope and get data
        # rm = visa.ResourceManager()
        # funcgen = rm.open_resource(NAME, timeout=20, chunk_size=1024000)
        #
        # # Get the frequency
        # freq = float(funcgen.query("SOUR1:FREQ?"))
        # print("returned frequency", freq)
        # funcgen.close()
        # rm.close()

        # query using the Signal Hound spectrum analyzer
        NUM_CAPTURES = 1
        DECIMATION = 1

        center, span, ref_level_dBm = [self.Rigol_carrier_frequency, 2.2*self.Rigol_FM_deviation, 20]

        # Open device
        handle = sa_open_device()["handle"]

        # Configure device
        sa_config_center_span(handle, center, span)
        sa_config_level(handle, ref_level_dBm)
        sa_config_IQ(handle, DECIMATION, 100.0e3);

        # Initialize
        sa_initiate(handle, SA_IQ, 0);
        return_len = sa_query_stream_info(handle)["return_len"]

        # Stream IQ
        print("Streaming...")
        sample_count = 0

        data = sa_get_IQ_32f(handle)["iq"]

        # Close device
        sa_close_device(handle)

        # get frequency of the peak
        iq = data

        # Create window
        window = get_window("hamming", len(iq))
        # Normalize window
        window *= len(window) / sum(window)
        # Window, FFT, normalize FFT output
        iq_data_FFT = numpy.fft.fftshift(numpy.fft.fft(iq * window) / len(window))
        spectrum = (iq_data_FFT.real ** 2 + iq_data_FFT.imag ** 2) / max(abs(iq_data_FFT) ** 2)

        freq_pts = np.linspace(center - span / 2, center + span / 2, len(iq_data_FFT))
        peak_freq = freq_pts[np.argmax(spectrum)]
        # print("resolution:", freq_pts[1] - freq_pts[0], "Hz")
        # print("peak:", peak_freq, "Hz")

        return peak_freq

    @kernel
    def set_Rigol_frequency_by_voltage(self, f_Rigol_modulation: TFloat):
        self.core.reset()

        # this fails because visa locks up the resource?

        # conversion to voltage assuming a linear response. this is not exactly true.
        # check the frequency on the oscilloscope. Rigol_modulation_volts is declared in BaseExperiment
        Rigol_modulation_volts = ((f_Rigol_modulation - self.Rigol_carrier_frequency)
                                       * 5 / self.Rigol_FM_deviation)
        self.zotino0.write_dac(4, Rigol_modulation_volts)
        self.zotino0.load()
        delay(0.1*ms)

    def run(self):

        for i,f in enumerate(self.scan_sequence1):
            print("request frequency",f)
            self.set_Rigol_frequency_by_voltage(f)
            self.actual_frequency_list[i] = self.query_Rigol_frequency()

        self.set_dataset("actual_frequencies", self.actual_frequency_list, archive=True)
        self.set_dataset("requested_frequencies", self.scan_sequence1, archive=True)
