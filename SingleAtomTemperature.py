"""
Two-shot experiment with a varying time between the shots over which the trap is turned off before being turned back
on just before the second shot
"""

from artiq.experiment import *
import csv
import numpy as np
from datetime import datetime as dt

from utilities.BaseExperiment import BaseExperiment
from subroutines.experiment_functions import atom_loading_2_experiment


class SingleAtomTemperature(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.setattr_argument("n_measurements", NumberValue(100, type='int', scale=1, ndecimals=0, step=1))
        # this is an argument for using a scan package, maybe
        self.scan_datasets = ["t_FORT_drop_sequence"]
        try:
            for dataset in self.scan_datasets:
                value = self.get_dataset(dataset)
                self.setattr_argument(dataset, StringValue(value))
        except KeyError as e:
            print(e)
            self.setattr_argument("t_FORT_drop_sequence", StringValue(
                'np.array([1.0, 10.0, 50.0, 100.])*us'))

        self.setattr_argument("no_first_shot", BooleanValue(False))
        self.setattr_argument("do_PGC_in_MOT", BooleanValue(False))
        self.setattr_argument("bins", NumberValue(50, ndecimals=0, step=1), "Histogram setup (set bins=0 for auto)")

        self.base.set_datasets_from_gui_args()
        print("build - done")

    def prepare(self):
        """
        performs initial calculations and sets parameter values before
        running the experiment. also sets data filename for now.

        any conversions from human-readable units to machine units (mu) are done here
        """
        self.base.prepare()

        self.t_exp_trigger = 1*ms

        # self.sampler_buffer = np.full(8, 0.0)
        # self.cooling_volts_ch = 7

        self.t_FORT_drop_list = eval(self.t_FORT_drop_sequence)
        self.n_iterations = len(self.t_FORT_drop_list)

        print("prepare - done")

    @kernel
    def run(self):
        self.base.initialize_hardware()
        self.base.initialize_datasets()
        self.expt()
        print("*************   Experiment finished   *************")

    @kernel
    def expt(self):
        """
        The experiment loop.

        :return:
        """
        self.zotino0.set_dac([0.0, 0.0, 0.0, 0.0],  # voltages must be floats or ARTIQ complains
                             channels=self.coil_channels)

        # todo: these are going to be regularly used, so put these in the base experiment
        self.set_dataset("SPCM0_RO1", [0], broadcast=True)
        self.set_dataset("SPCM0_RO2", [0], broadcast=True)
        self.set_dataset("photocount_bins", [self.bins], broadcast=True)


        # turn on cooling MOT AOMs
        self.dds_cooling_DP.sw.on() # cooling double pass
        self.dds_AOM_A2.sw.on()
        self.dds_AOM_A3.sw.on()
        self.dds_AOM_A1.sw.on()
        self.dds_AOM_A6.sw.on()
        self.dds_AOM_A4.sw.on()
        self.dds_AOM_A5.sw.on()
        delay(1 * ms)

        # delay(2000*ms) # wait for AOMS to thermalize in case they have been off.

        if self.enable_laser_feedback:
            self.laser_stabilizer.run()
        delay(1*ms)

        SPCM0_RO1 = 0
        SPCM0_RO2 = 0

        iteration = 0
        for t_FORT_drop in self.t_FORT_drop_list:

            ### These are the datasets for plotting only. We restart them each iteration
            self.set_dataset("SPCM0_RO1_current_iteration", [0], broadcast=True)
            self.set_dataset("SPCM0_RO2_current_iteration", [0], broadcast=True)

            self.set_dataset("SPCM1_RO1_current_iteration", [0], broadcast=True)
            self.set_dataset("SPCM1_RO2_current_iteration", [0], broadcast=True)

            self.t_FORT_drop = t_FORT_drop
            atom_loading_2_experiment(self)

            # ### loop the experiment sequence
            # for measurement in range(self.n_measurements):
            #
            #     if self.enable_laser_feedback:
            #         # if measurement % 10 == 0:
            #         self.laser_stabilizer.run()
            #         delay(1 * ms)
            #         self.dds_FORT.sw.on()
            #         self.dds_FORT.set(frequency=self.f_FORT - 30 * MHz, amplitude=self.stabilizer_FORT.amplitude)
            #
            #     load_MOT_and_FORT(self)
            #
            #     # set the cooling DP AOM to the readout settings
            #     self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO, amplitude=self.ampl_cooling_DP_MOT)
            #
            #     if not self.no_first_shot:
            #         # take the first shot
            #         self.dds_cooling_DP.sw.on()
            #         t_gate_end = self.ttl0.gate_rising(self.t_SPCM_first_shot)
            #         SPCM0_RO1 = self.ttl0.count(t_gate_end)
            #         delay(1*ms)
            #         self.dds_cooling_DP.sw.off()
            #
            #     # turn the FORT off
            #     self.dds_FORT.set(frequency=self.f_FORT - 30 * MHz, amplitude=self.stabilizer_FORT.amplitude)
            #
            #     delay(t_FORT_drop)
            #
            #     # turn the FORT on
            #     self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude)
            #     delay(1*ms)
            #
            #     # take the second shot
            #     self.dds_cooling_DP.sw.on()
            #     t_gate_end = self.ttl0.gate_rising(self.t_SPCM_second_shot)
            #     SPCM0_RO2 = self.ttl0.count(t_gate_end)
            #     delay(1*ms)
            #     self.dds_cooling_DP.sw.off()
            #
            #     self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
            #
            #     delay(2*ms)
            #
            #     iteration += 1
            #
            #     # update the datasets
            #     if not self.no_first_shot:
            #         self.append_to_dataset('SPCM0_RO1', SPCM0_RO1)
            #         self.append_to_dataset('SPCM0_RO1_current_iteration', SPCM0_RO1)
            #
            #     # update the datasets
            #     self.append_to_dataset('SPCM0_RO2', SPCM0_RO2)
            #     self.append_to_dataset('SPCM0_RO2_current_iteration', SPCM0_RO2)
            #     self.set_dataset("iteration", iteration, broadcast=True)

        delay(1*ms)

        ### leave MOT on at end of experiment, but turn off the FORT
        self.dds_cooling_DP.sw.on()
        self.zotino0.set_dac([self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
                             channels=self.coil_channels)
        ### effectively turn the FORT AOM off
        self.dds_FORT.set(frequency=self.f_FORT - 30 * MHz, amplitude=self.stabilizer_FORT.amplitude)
        # set the cooling DP AOM to the MOT settings

        ### finally, in case the worker refuses to die
        self.write_results()