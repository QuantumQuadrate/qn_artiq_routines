"""
A simple experiment to scan a variable. This avoids generalVariableScan and is a simple code.
In the example below, I am scanning t_FORT_drop time which can be used for atom temperature measurement. For this,
I set self.t_FORT_drop = scan_variable in the loop and run atom_loading_2_experiment:

"""

from artiq.experiment import *
import numpy as np
from datetime import datetime as dt

import sys, os
### get the current working directory
current_working_directory = os.getcwd()
cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")

from utilities.BaseExperiment import BaseExperiment
from subroutines.experiment_functions import atom_loading_2_experiment

class Scan_a_variable(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.setattr_argument("n_measurements", NumberValue(100, type='int', scale=1, ndecimals=0, step=1))

        ### this is an argument for using a scan package, maybe
        self.scan_datasets = ["scan_sequence"]
        try:
            for dataset in self.scan_datasets:
                value = self.get_dataset(dataset)
                self.setattr_argument(dataset, StringValue(value))
        except KeyError as e:
            print(e)
            self.setattr_argument("scan_sequence", StringValue(
                'np.array([1.0, 10.0, 100.])*us'))

        self.setattr_argument("no_first_shot", BooleanValue(False))

        self.base.set_datasets_from_gui_args()
        print("build - done")

    def prepare(self):
        """
        performs initial calculations and sets parameter values before
        running the experiment. also sets data filename for now.

        any conversions from human-readable units to machine units (mu) are done here
        """
        self.base.prepare()

        self.scan_list = eval(self.scan_sequence)
        self.n_iterations = len(self.scan_list)

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

        ### Set the coils to MOT loading setting
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
            channels=self.coil_channels)

        ### set the cooling DP AOM to the MOT settings
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)

        self.dds_cooling_DP.sw.on()  ### turn on cooling
        self.ttl_repump_switch.off()  ### turn on MOT RP

        self.dds_AOM_A1.sw.on()
        self.dds_AOM_A2.sw.on()
        self.dds_AOM_A3.sw.on()
        self.dds_AOM_A4.sw.on()
        delay(0.1 * ms)
        self.dds_AOM_A5.sw.on()
        self.dds_AOM_A6.sw.on()

        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude)
        self.dds_FORT.sw.on()

        delay(1 * ms)


        if self.enable_laser_feedback:
            self.laser_stabilizer.run()
        delay(1*ms)


        for scan_variable in self.scan_list:
            self.t_FORT_drop = scan_variable

            ### These are the datasets for plotting only. but necessary if use end_measurement. We restart them each iteration
            self.set_dataset("SPCM0_RO1_current_iteration", [0], broadcast=True)
            self.set_dataset("SPCM0_RO2_current_iteration", [0], broadcast=True)

            self.set_dataset("SPCM1_RO1_current_iteration", [0], broadcast=True)
            self.set_dataset("SPCM1_RO2_current_iteration", [0], broadcast=True)

            atom_loading_2_experiment(self)


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