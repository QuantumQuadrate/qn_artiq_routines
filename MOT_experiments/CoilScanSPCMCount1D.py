"""
This code scans the coils to find the optimum coil parameters that put the MOT at the focus of the
parabolic mirror. It saves the photon counts in a file.

"""
import os

from artiq.experiment import *
import csv
from datetime import datetime as dt
import numpy as np
<<<<<<< HEAD
import sys, os
# get the current working directory
current_working_directory = os.getcwd()
cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")
=======
import sys
sys.path.append('C:\\Users\\jakeuribe\\artiq-master\\repository\\qn_artiq_routines\\')
sys.path.append('C:\\Users\\jakeuribe\\artiq-master\\repository\\')
>>>>>>> 989bec7 (unshelved)

from utilities.BaseExperiment import BaseExperiment



class CoilScanSPCMCount1D(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.setattr_argument("print_time_estimate_only", BooleanValue(True))

        self.scan_datasets = ["coil_V_array"]
        group = "Coil steps"

        self.setattr_argument("coil_name", EnumerationValue(['AZ_bottom','AZ_top','AX','AY']))

        try:
            for dataset in self.scan_datasets:
                value = self.get_dataset(dataset)
                self.setattr_argument(dataset, StringValue(value), group)
        except KeyError as e:
            print(e)
            self.setattr_argument("coil_V_array", StringValue(
                '[6 - l*(6 - 8)/20 for l in range(20)]'))

        self.setattr_argument("differential_scan", BooleanValue(True))

        self.setattr_argument("iterations", NumberValue(10, ndecimals=0, step=1, scale=1, type='int'))

        # for now, just run laser feedback at beginning
        self.setattr_argument("enable_laser_feedback", BooleanValue(True), "Laser feedback")

        # group1 = "Coil tune settings"
        # self.setattr_argument("enable_coil_tune", BooleanValue(True), group1)
        #
        # self.setattr_argument("Acceptable_SPCM_Count", NumberValue(10000.0, ndecimals=1, step=1), group1)
        # self.setattr_argument("FORT_AOM_on", BooleanValue(False), group1)
        # self.setattr_argument("coil_volts_multiplier",
        #                       NumberValue(3.3), group1)  # scales the value read by the Sampler
        # self.setattr_argument("differential_mode",
        #                       BooleanValue(False), group1)  # scan the coils with respect to the current settings
        # self.setattr_argument("differential_multiplier",
        #                       NumberValue(0.4), group1)  # scales the value read by the Sampler

        # # exposure time of the SPCM
        # self.setattr_argument("Tune_SPCM_exposure", NumberValue(300 * ms, unit='ms'), group1)
        # # saturation limit of the SPCM in counts/s. Can be increased to 10**7 safely, but not higher than 3*10**7.
        # self.setattr_argument("sat1s", NumberValue(1 * 10 ** 5), group1)  # saturation limit in counts/dt.
        # self.setattr_argument("print_count_rate", BooleanValue(True), group1)

        group2 = "Data options"
        self.setattr_argument("datafile", StringValue('coil_scan_1D.csv'), group2)
        self.setattr_argument("prepend_date_to_datafile", BooleanValue(True), group2)
        self.setattr_argument("print_measurement_number", BooleanValue(False), group2)
        self.setattr_argument("print_meas_result", BooleanValue(False), group2)
        self.setattr_argument("save_data", BooleanValue(True), group2)

        self.base.set_datasets_from_gui_args()
        print("build - done")

    def prepare(self):
        """
        performs initial calculations and sets parameter values before
        running the experiment. also sets data filename for now.

        any conversions from human-readable units to machine units (mu) are done here
        """

        self.base.prepare()


        # where to store the data
        self.t_experiment_run = dt.now().strftime("%Y%m%d_%H%M%S")
        self.datadir = os.getcwd()
        if self.prepend_date_to_datafile:
            self.datafile = os.path.join(self.datadir, self.t_experiment_run + '_' + self.datafile)
        else:
            self.datafile = os.path.join(self.datadir + self.datafile)

        # save a copy of the string defining the scan variable
        self.V_array_string = self.coil_V_array

        # the default coil volts for the channel we're going to change
        if self.coil_name == 'AZ_bottom':
            self.coil_index = 0
            # self.coil_V0 = self.AZ_bottom_volts_MOT
        elif self.coil_name == 'AZ_top':
            self.coil_index = 1
            # self.coil_V0 = self.AZ_bottom_volts_MOT
        elif self.coil_name == 'AX':
            self.coil_index = 2
            # self.coil_V0 = self.AX
        else:
            self.coil_index = 3
            # self.coil_V0 = self.AY

        # evaluate the strings we used to define the coil steps in the GUI.
        self.coil_V_array = eval(self.coil_V_array)
        self.V_steps = len(self.coil_V_array)
        self.default_volts = [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT]
        self.coil_volts = [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT]

        if self.differential_scan:
            self.coil_V_array += self.coil_volts[self.coil_index]

        self.sampler_buffer = np.full(8, 0.0)
        self.control_volts_channels = [0, 1, 2, 3]  # the sampler channels to read

        self.count_rate_dataset = 'photocounts_per_s'
        self.set_dataset(self.count_rate_dataset,
                         [0.0],
                         broadcast=True)

        self.maxcount_dataset = 'SPCM_max_counts_and_volts'
        self.set_dataset(self.count_rate_dataset,
                         [0.0, *self.default_volts],
                         broadcast=True)

        print("prepare - done")

    # @rpc(flags={'async'})
    # def update_coil_values(self)->TArray(TFloat):
    #     self.AZ_bottom_volts_MOT = self.get_dataset("AZ_bottom_volts_MOT")
    #     self.AZ_top_volts_MOT = self.get_dataset("AZ_top_volts_MOT")
    #     self.AX_volts_MOT = self.get_dataset("AX_volts_MOT")
    #     self.AY_volts_MOT = self.get_dataset("AY_volts_MOT")
    #     self.coil_volts = [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT]
    #     self.coil_V0 = [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT][
    #         self.coil_index]
    #     # if self.differential_scan:
    #     #     self.coil_V_array += self.coil_volts[self.coil_index]
    #     # print(self.coil_V_array)
    #     # assert max(abs(self.coil_V_array)) <= 10, "Zotino can not output > 10V. Decrease V scan range"
    #     coil_V_array = self.coil_V_array + self.coil_volts[self.coil_index]
    #     return coil_V_array

    @kernel
    def run(self):
        self.base.initialize_hardware()
        # if self.enable_coil_tune:
        #     self.tune_experiment()
        self.scan_experiment()
        print("Experiment finished.")

    # @kernel
    # def tune_experiment(self):
    #
    #     # self.update_coil_values()
    #
    #     ### Allows you to scan the coils to maximize the SPCM counts before running the atom loading exp (LoadExpt).
    #     # Turn on AOMs to load the MOT.
    #     self.dds_cooling_DP.sw.on()
    #
    #
    #     self.dds_AOM_A2.sw.on()
    #     self.dds_AOM_A3.sw.on()
    #     self.dds_AOM_A1.sw.on()
    #     self.dds_AOM_A6.sw.on()
    #     self.dds_AOM_A4.sw.on()
    #     self.dds_AOM_A5.sw.on()
    #
    #     # wait for AOMs to thermalize
    #     delay(3000 * ms)
    #
    #     if self.FORT_AOM_on:
    #         self.dds_FORT.sw.on()
    #     else:
    #         self.dds_FORT.sw.off()
    #
    #     control_volts = [0.0] * 4
    #     best_volts = [0.0] * 4  # the volts that we had when the best MOT was found
    #     max_count_rate = 0.0
    #     count_rate_Hz = 1.0
    #
    #     delay(10 * ms)
    #
    #     # for outputting a voltage proportional to counts read from the SPCM on ttl0
    #     ch = 6
    #     self.zotino0.write_dac(ch, 0.0)
    #     self.zotino0.load()
    #
    #     Satdt = self.sat1s * self.Tune_SPCM_exposure  # saturation limit in counts/dt.
    #     delay(1 * ms)
    #
    #     if self.enable_laser_feedback:
    #         self.laser_stabilizer.run()
    #         if self.FORT_AOM_on:
    #             self.dds_FORT.sw.on()
    #
    #     print("****************   ready for tuning   ********************")
    #
    #     i = 0
    #     while count_rate_Hz < self.Acceptable_SPCM_Count:
    #
    #         if self.enable_laser_feedback:
    #             if (i % 100) == 0:
    #                 print("running feedback")
    #                 self.core.break_realtime()
    #                 self.laser_stabilizer.run()
    #                 delay(100 * ms)
    #                 if self.FORT_AOM_on:
    #                     self.dds_FORT.sw.on()
    #
    #         i = i + 1
    #
    #         tend1 = self.ttl0.gate_rising(self.Tune_SPCM_exposure)
    #         count1 = self.ttl0.count(tend1)
    #         count_rate_Hz = count1 / self.Tune_SPCM_exposure
    #         if self.print_count_rate:
    #             print(round(count_rate_Hz))
    #         delay(10 * ms)
    #
    #         volt1 = count1 * 5 / Satdt  # the voltage from zotino0, port 7. Saturation limit corresponds to 5V.
    #         self.zotino0.write_dac(ch, volt1)
    #         self.zotino0.load()
    #         delay(2 * ms)
    #
    #         self.append_to_dataset(self.count_rate_dataset, count_rate_Hz)
    #         if count_rate_Hz > max_count_rate:
    #             max_count_rate = count_rate_Hz
    #             best_volts = control_volts
    #             self.set_dataset(self.maxcount_dataset, [count_rate_Hz] + best_volts)
    #         delay(1 * ms)
    #
    #         self.sampler1.sample(self.sampler_buffer)
    #         if self.differential_mode:
    #             control_volts = [self.sampler_buffer[ch] * self.differential_multiplier + self.default_volts[ch]
    #                              for ch in self.control_volts_channels]
    #         else:
    #             control_volts = [self.sampler_buffer[ch] * self.coil_volts_multiplier
    #                              for ch in self.control_volts_channels]
    #
    #         delay(1 * ms)
    #
    #         # set coils based on the sampler values we read
    #         self.zotino0.set_dac(
    #             control_volts,
    #             channels=self.coil_channels)
    #
    #     volt_datasets = ["AZ_bottom_volts_MOT", "AZ_top_volts_MOT", "AX_volts_MOT", "AY_volts_MOT"]
    #     for i in range(4):
    #         self.set_dataset(volt_datasets[i], best_volts[i], broadcast=True)
    #
    #     print("Best volts [VZ_bottom,VZ_top,Vx,Vy]:")
    #     print(best_volts)
    #     print("Best count rate (kHz):")
    #     print(max_count_rate / 1e3)
    #
    #     delay(100 * ms)

    @kernel
    def scan_experiment(self):

        print("estimated duration of scan:",
              self.iterations * self.V_steps * (self.t_MOT_loading + self.t_SPCM_exposure) / 60,
              "minutes")
        if not self.print_time_estimate_only:

            # self.update_coil_values()

            self.file_setup(rowheaders=['counts','AZ_bottom V','AZ_top V','AY V','AX V','cooling PD V'])

            delay(10 * ms)
            # the starting point
            self.zotino0.set_dac([self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
                                 channels=self.coil_channels)
            delay(1 * ms)  # avoid RTIOSequence error

            # Turn on AOMs to load the MOT.
            self.dds_cooling_DP.sw.on()


            self.dds_AOM_A2.sw.on()
            self.dds_AOM_A3.sw.on()
            self.dds_AOM_A1.sw.on()
            self.dds_AOM_A6.sw.on()
            self.dds_AOM_A4.sw.on()
            self.dds_AOM_A5.sw.on()

            # wait for AOMs to thermalize
            delay(4000 * ms)

            step = 0

            delay(100*ms)
            if self.enable_laser_feedback:
                self.laser_stabilizer.run()
                delay(100 * ms)

            # start with the coil defaults as a check
            self.zotino0.set_dac(
                self.coil_volts,
                channels=self.coil_channels)

            delay(1000*ms)

            for i in range(self.iterations):
                if self.enable_laser_feedback:
                    print("running feedback")
                    self.core.break_realtime()
                    self.laser_stabilizer.run()
                    delay(100 * ms)

                for V in self.coil_V_array:

                    # do the experiment sequence
                    delay(10 * ms)

                    # update coil values
                    self.coil_volts[self.coil_index] = V

                    delay(1*ms)
                    self.zotino0.set_dac(
                        self.coil_volts,
                        channels=self.coil_channels)

                    # wait for the MOT to load
                    delay_mu(self.t_MOT_loading_mu)

                    # take the shot
                    t_gate_end = self.ttl0.gate_rising(self.t_SPCM_exposure)
                    counts = self.ttl0.count(t_gate_end)

                    ## trigger the Luca
                    delay(1*ms)
                    self.ttl6.pulse(5*ms)
                    delay(60*ms)

                    if self.print_meas_result:
                        print("counts", counts)
                    delay(1 * ms)
                    if self.save_data:
                        ### all items in the list must be the same type that is why 0.0 is added to counts
                        self.sampler0.sample(self.sampler_buffer) # check cooling laser power
                        self.file_write([counts+0.0,
                                         self.coil_volts[0],
                                         self.coil_volts[1],
                                         self.coil_volts[2],
                                         self.coil_volts[3],
                                         self.sampler_buffer[self.cooling_volts_ch]])
                    delay(10 * ms)
                    step += 1

            ### reset parameters
            delay(10*ms)

            self.zotino0.set_dac([0.0, 0.0, 0.0, 0.0],
                                 channels=self.coil_channels)

            self.cleanup()  # asynchronous stuff like closing files


    @rpc(flags={"async"})  # means this code runs asynchronously; won't block the rtio counter
    def cleanup(self):
        self.file_obj.close()

    @rpc(flags={"async"})
    def file_setup(self, rowheaders=[]):
        # open file once, then close it at end of the experiment
        self.file_obj = open(self.datafile, 'a', newline='')
        self.csvwriter = csv.writer(self.file_obj)
        if rowheaders != []:
            self.csvwriter.writerow(rowheaders)

    @rpc(flags={"async"})
    def file_write(self, data):
        self.csvwriter.writerow(data)

