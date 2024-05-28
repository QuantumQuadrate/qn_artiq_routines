from artiq.experiment import *

import numpy as np
import os
import cv2
from PIL import Image
from thorlabs_tsi_sdk.tl_camera import TLCameraSDK, OPERATION_MODE
import matplotlib.pyplot as plt
from datetime import datetime as dt

import sys
# get the current working directory
current_working_directory = os.getcwd()
cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")
from utilities.BaseExperiment import BaseExperiment

# this is where your experiment function should live
from subroutines.experiment_functions import *
import subroutines.experiment_functions as exp_functions

class ExperimentCycler(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        # the number of measurements to be made for a certain setting of the
        # experiment parameters
        self.setattr_argument("n_measurements", NumberValue(10, ndecimals=0, step=1))

        experiment_function_names_list = [x for x in dir(exp_functions)
            if ('__' not in x and str(type(getattr(exp_functions,x)))=="<class 'function'>"
                and 'experiment' in x)]

        # a function that take no arguments that gets imported and run
        # self.setattr_argument('experiment_function', StringValue('test'))
        self.setattr_argument('experiment_function',
                              EnumerationValue(experiment_function_names_list))

        self.base.set_datasets_from_gui_args()
        print("build - done")

    def prepare(self):
        """
        performs initial calculations and sets parameter values before
        running the experiment.
        """
        self.base.prepare()

        try:
            self.experiment_name = self.experiment_function
            self.experiment_function = lambda :eval(self.experiment_name)(self)
        except NameError as e:
            print(f"The function {self.experiment_name} is not defined. Did you forget to import it?")
            raise

        self.measurement = 0
        self.counts = 0
        self.counts2 = 0

    @kernel
    def hardware_init(self):
        self.base.initialize_hardware()

    # todo: this should really be determined by the specific experiment eventually
    def initialize_datasets(self):
        self.set_dataset("n_measurements", self.n_measurements, broadcast=True)
        self.set_dataset("photocounts", [0], broadcast=True)
        self.set_dataset("photocounts2", [0], broadcast=True)
        self.set_dataset("photocount_bins", [50], broadcast=True)
        self.set_dataset("iteration", 0, broadcast=True)

        self.set_dataset("feedbackchannels",
                         [ch.name for ch in self.laser_stabilizer.all_channels],
                         broadcast=True, persist=True)

        value = 0.0
        for ch_i in range(len(self.laser_stabilizer.all_channels)):
            self.set_dataset(self.laser_stabilizer.all_channels[ch_i].dB_history_dataset,
                             [float(self.initial_RF_dB_values[ch_i])], broadcast=True, persist=True)

    def reset_datasets(self):
        """
        set datasets that are redefined each iteration.
        :return:
        """
        self.set_dataset("test_dataset", [0], broadcast=True)

        # typically these datasets are used for plotting which would be meaningless if we continued to append to the data,
        # e.g. for the second readout histogram which we expect in general will change as experiment parameters induce
        # different amount of atom loss.
        self.set_dataset('photocounts_current_iteration', [0], broadcast=True)
        self.set_dataset('photocounts2_current_iteration', [0], broadcast=True)

        # no reason to let these datasets grow to huge lengths
        for ch in self.laser_stabilizer.all_channels:
            self.set_dataset(ch.dataset, [self.get_dataset(ch.dataset)[-1]], broadcast=True)

    def run(self):
        """
        Loop the experiment function, pausing only for higher-priority experiments
        in the same pipeline
        """

        self.initialize_datasets()

        iteration = 0

        while True:

            self.hardware_init()
            self.reset_datasets()

            # the measurement loop.
            self.experiment_function()
            self.write_results({'name': self.experiment_name[:-11]})

            if self.scheduler.check_pause():
                self.core.comm.close()  # put the hardware in a safe state before checking pause
                self.scheduler.pause()  # check if we need to run a new experiment*

                # todo: build executes but doesn't seem to update the parameters used by the
                #  experiment. I tried updating n_measurements
                # after pause is done, we want to re-initialize variables in case they have changed
                # self.base.build()

            iteration += 1
            self.set_dataset("iteration", iteration, broadcast=True)

            # todo: add in compute loading so we can log the rate and retention





