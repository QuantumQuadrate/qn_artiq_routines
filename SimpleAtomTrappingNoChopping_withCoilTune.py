"""
A simple experiment to look for a trapped single atom signal

1. Turn on cooling and RP AOMs
Experiment cycle (repeats n times)
2. Turn on magnetic fields (Zotino)
3. Turn on fiber AOMs (Urukul), wait some time to load the MOT
4. Turn on the dipole trap AOM
5. Turn off the quadrupole fields (PGC phase)
6. Turn off the fiber AOMs
7. Turn on the fiber AOMs and read from a single SPCM (TTL) for a certain exposure time
8. Store the number of counts registered by the SPCM in an array
End of experiment
9. Save the array of counts to a file


* Tested that all the delays and timings are functioning as expected by pulsing TTLs at different points and monitoring
    with an oscilloscope.

"""

from artiq.experiment import *
import csv
import numpy as np
from datetime import datetime as dt

from utilities.BaseExperiment import BaseExperiment


class SimpleAtomTrapNoChopCoilTune(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        group1 = "Coil tune settings"
        self.setattr_argument("Acceptable_SPCM_Count", NumberValue(10000.0, ndecimals=1, step=1), group1)
        self.setattr_argument("FORT_AOM_on", BooleanValue(False), group1)
        self.setattr_argument("coil_volts_multiplier",
                              NumberValue(3.3), group1)  # scales the value read by the Sampler
        self.setattr_argument("differential_mode",
                              BooleanValue(False), group1)  # scan the coils with respect to the current settings
        self.setattr_argument("differential_multiplier",
                              NumberValue(0.4), group1)  # scales the value read by the Sampler

        # exposure time of the SPCM
        self.setattr_argument("Tune_SPCM_exposure", NumberValue(300 * ms, unit='ms'), group1)
        # saturation limit of the SPCM in counts/s. Can be increased to 10**7 safely, but not higher than 3*10**7.
        self.setattr_argument("sat1s", NumberValue(1 * 10 ** 5), group1)  # saturation limit in counts/dt.
        self.setattr_argument("print_count_rate", BooleanValue(True), group1)


        group2 = "Atom loading settings"
        self.setattr_argument("n_measurements", NumberValue(10, ndecimals=0, step=1), group2)
        self.setattr_argument("print_measurement_number", BooleanValue(False), group2)
        self.setattr_argument("Load_SPCM_exposure", NumberValue(10 * ms, unit='ms'), group2)
        self.setattr_argument("expose_with_MOT_on", BooleanValue(False), group2)
        self.setattr_argument("bins", NumberValue(50, ndecimals=0, step=1), group2)
        self.setattr_argument("print_counts", BooleanValue(True), group2)
        self.setattr_argument("control_experiment", BooleanValue(False), group2)

        # when to run the AOM feedback (after how many iterations in the for loops)
        group3 = "MOT beams feedback"
        self.setattr_argument("enable_laser_feedback", BooleanValue(default=True), group3)
        self.setattr_argument("AOM_feedback_period_cycles", NumberValue(200), group3)

        self.base.set_datasets_from_gui_args()
        print("build - done")

    def prepare(self):
        """
        performs initial calculations and sets parameter values before
        running the experiment. also sets data filename for now.

        any conversions from human-readable units to machine units (mu) are done here
        """
        self.base.prepare()

        self.sampler_buffer = np.zeros(8)

        self.control_volts_channels = [0, 1, 2, 3]  # the sampler channels to read
        self.default_volts = [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT]

        self.count_rate_dataset = 'SPCM_count_rate_Hz'
        self.set_dataset(self.count_rate_dataset,
                         [0.0],
                         broadcast=True)

        self.maxcount_dataset = 'SPCM_max_counts_and_volts'
        self.set_dataset(self.count_rate_dataset,
                         [0.0, *self.default_volts],
                         broadcast=True)

        self.t_exp_trigger = 1*ms

        self.sampler_buffer = np.full(8, 0.0)
        self.cooling_volts_ch = 7

        # self.hist_bins = np.zeros(self.bins, dtype=int)
        # self.photocounts = np.full(self.n_measurements, 0.0)

        print("prepare - done")

    @kernel
    def run(self):
        self.base.initialize_hardware()
        self.TuneExpt()
        self.LoadExpt()
        print("Experiment finished.")

    @kernel
    def TuneExpt(self):
        ### Allows you to scan the coils to maximize the SPCM counts before running the atom loading exp (LoadExpt).
        # Turn on AOMs to load the MOT.
        self.dds_cooling_DP.sw.on()
        self.dds_cooling_SP.sw.on()
        self.dds_MOT_RP.sw.on()
        self.dds_AOM_A2.sw.on()
        self.dds_AOM_A3.sw.on()
        self.dds_AOM_A1.sw.on()
        self.dds_AOM_A6.sw.on()
        self.dds_AOM_A4.sw.on()
        self.dds_AOM_A5.sw.on()

        # wait for AOMs to thermalize
        delay(3000 * ms)

        if self.FORT_AOM_on:
            self.dds_FORT.sw.on()
        else:
            self.dds_FORT.sw.off()


        control_volts = [0.0]*4
        best_volts = [0.0]*4 # the volts that we had when the best MOT was found
        max_count_rate = 0.0
        count_rate_Hz = 1.0

        delay(10 * ms)

        # for outputting a voltage proportional to counts read from the SPCM on ttl0
        ch = 6
        self.zotino0.write_dac(ch, 0.0)
        self.zotino0.load()

        Satdt = self.sat1s * self.Tune_SPCM_exposure  # saturation limit in counts/dt.
        delay(1 * ms)

        if self.enable_laser_feedback:
            self.laser_stabilizer.run()
            if self.FORT_AOM_on:
                self.dds_FORT.sw.on()


        print("****************   ready for tuning   ********************")

        # for i in range(self.n_steps):
        i = 0
        while count_rate_Hz < self.Acceptable_SPCM_Count:

            if self.enable_laser_feedback:
                if (i % self.AOM_feedback_period_cycles) == 0:
                    print("running feedback")
                    self.core.break_realtime()
                    self.laser_stabilizer.run()
                    delay(100*ms)
                    if self.FORT_AOM_on:
                        self.dds_FORT.sw.on()

            i = i + 1

            tend1 = self.ttl0.gate_rising(self.Tune_SPCM_exposure)
            count1 = self.ttl0.count(tend1)
            count_rate_Hz = count1 / self.Tune_SPCM_exposure
            if self.print_count_rate:
                print(round(count_rate_Hz))
            delay(10 * ms)

            volt1 = count1 * 5 / Satdt  # the voltage from zotino0, port 7. Saturation limit corresponds to 5V.
            self.zotino0.write_dac(ch, volt1)
            self.zotino0.load()
            delay(2 * ms)


            self.append_to_dataset(self.count_rate_dataset, count_rate_Hz)
            if count_rate_Hz > max_count_rate:
                max_count_rate = count_rate_Hz
                best_volts = control_volts
                self.set_dataset(self.maxcount_dataset,[count_rate_Hz]+best_volts)
            delay(1 * ms)

            self.sampler1.sample(self.sampler_buffer)
            if self.differential_mode:
                control_volts = [self.sampler_buffer[ch]*self.differential_multiplier + self.default_volts[ch]
                                 for ch in self.control_volts_channels]
            else:
                control_volts = [self.sampler_buffer[ch] * self.coil_volts_multiplier
                                 for ch in self.control_volts_channels]

            delay(1*ms)

            # set coils based on the sampler values we read
            self.zotino0.set_dac(
                control_volts,
                channels=self.coil_channels)


        volt_datasets = ["AZ_bottom_volts_MOT", "AZ_top_volts_MOT", "AX_volts_MOT", "AY_volts_MOT"]
        for i in range(4):
            self.set_dataset(volt_datasets[i], best_volts[i], broadcast=True)

        print("Best volts [VZ_bottom,VZ_top,Vx,Vy]:")
        print(best_volts)
        print("Best count rate (kHz):")
        print(max_count_rate/1e3)

        delay(100*ms)

    @rpc(flags={'async'})
    def update_coil_values(self):
        self.AZ_bottom_volts_MOT = self.get_dataset("AZ_bottom_volts_MOT")
        self.AZ_top_volts_MOT = self.get_dataset("AZ_top_volts_MOT")
        self.AX_volts_MOT = self.get_dataset("AX_volts_MOT")
        self.AY_volts_MOT = self.get_dataset("AY_volts_MOT")

    @kernel
    def LoadExpt(self):
        """
        The experiment loop.

        :return:
        """

        self.update_coil_values()

        self.set_dataset("photocounts", [0], broadcast=True)
        self.set_dataset("photocount_bins", [self.bins], broadcast=True)

        # loop the experiment sequence
        for measurement in range(self.n_measurements):

            if self.enable_laser_feedback:
                if (measurement % self.AOM_feedback_period_cycles) == 0:
                    print("running feedback")
                    self.core.break_realtime()
                    self.laser_stabilizer.run()
                    delay(100 * ms)

            self.ttl7.pulse(self.t_exp_trigger) # in case we want to look at signals on an oscilloscope

            if self.control_experiment and measurement % 2 == 0:
                # change the Y MOT coil enough to lose the MOT but not enough to significantly change the fluorescence
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT-0.8],
                    channels=self.coil_channels)
            else:
                # Set magnetic fields for MOT loading
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
                    channels=self.coil_channels)
            delay(1 * ms)  # avoid RTIOSequence error

            # wait for the MOT to load
            delay_mu(self.t_MOT_loading_mu)

            # turn on the dipole trap and wait to load atoms
            self.dds_FORT.sw.on()
            delay_mu(self.t_FORT_loading_mu)

            if not self.expose_with_MOT_on:
                # change the magnetic fields for imaging
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_RO, self.AZ_top_volts_RO, self.AX_volts_RO, self.AY_volts_RO],
                    channels=self.coil_channels)
                delay(10 * ms)

                # change AOMs to "imaging" settings
                # self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.ampl_FORT_RO)
                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO, amplitude=self.ampl_cooling_DP_RO)
                delay(1000 * ms)

            # take the shot
            t_gate_end = self.ttl0.gate_rising(self.Load_SPCM_exposure)
            counts = self.ttl0.count(t_gate_end)/self.Load_SPCM_exposure
            delay(1*ms)
            if self.print_counts:
                print(counts)
            delay(10 * ms)


            self.append_to_dataset('photocounts', counts)

            if self.print_measurement_number:
                print("measurement", measurement)
            delay(10*ms)

            self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)

        delay(1*ms)
        # leave MOT on at end of experiment, but turn off the FORT
        self.dds_FORT.sw.off()

        self.zotino0.set_dac([self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
                             channels=self.coil_channels)