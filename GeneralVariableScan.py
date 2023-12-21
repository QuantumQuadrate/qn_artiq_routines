from artiq.experiment import *

import numpy as np
import os
import cv2
from PIL import Image
from thorlabs_tsi_sdk.tl_camera import TLCameraSDK, OPERATION_MODE
import matplotlib.pyplot as plt
from datetime import datetime as dt

import sys
sys.path.append('C:\\Networking Experiment\\artiq codes\\artiq-master\\repository\\qn_artiq_routines\\')
from utilities.BaseExperiment import BaseExperiment
from subroutines.experiment_functions import *

class GeneralVariableScan(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        # the number of measurements to be made for a certain setting of the
        # experiment parameters
        self.setattr_argument("n_measurements", NumberValue(100, ndecimals=0, step=1))
        self.setattr_argument('scan_variable1', StringValue('t_blowaway'))
        # todo: need some error handling, or make this a drop box
        self.setattr_argument("scan_sequence1", StringValue(
            'np.array([0.000,0.005,0.02,0.05])*ms'))

        # this variable is optional
        self.setattr_argument('scan_variable2', StringValue(''))
        # todo: need some error handling, or make this a drop box
        self.setattr_argument("scan_sequence2", StringValue(
            'np.linspace(-2,2,5)*V'))

        # a function that take no arguments that gets imported and run
        self.setattr_argument('experiment_function', StringValue('test'))

        # toggles an interleaved control experiment, but what this means or whether
        # it has an effect depends on the experiment script
        self.setattr_argument("control_experiment", BooleanValue(False), "Control experiment")

        self.base.set_datasets_from_gui_args()
        print("build - done")

    def prepare(self):
        """
        performs initial calculations and sets parameter values before
        running the experiment.
        """
        self.base.prepare()

        self.scan_variable1 = str(self.scan_variable1)
        self.scan_variable2 = str(self.scan_variable2)

        assert hasattr(self,self.scan_variable1), (f"There is no ExperimentVariable "+self.scan_variable1+
                                                  ". Did you mistype it?")
        # todo: add some error handling
        # self.scan_variable1 = eval(self.scan_variable1) # probably neither necessary nor desired
        self.scan_sequence1 = eval(self.scan_sequence1)
        self.n_iterations1 = len(self.scan_sequence1) # this might not get used

        if self.scan_variable2 != '':
            assert hasattr(self, self.scan_variable2), (f"There is no ExperimentVariable " + self.scan_variable2 +
                                                        ". Did you mistype it?")
            # todo: add some error handling
            # self.scan_variable2 = eval(self.scan_variable2)
            self.n_iterations2 = len(self.scan_sequence2)
        else:
            self.scan_sequence2 = np.zeros(1) # we need a non-empty array to have a definite type for ARTIQ
            self.n_iterations2 = 0

        try:
            experiment_name = self.experiment_function
            self.experiment_function = lambda :eval(experiment_name)(self)
        except NameError as e:
            print(f"The function {experiment_name} is not defined. Did you forget to import it?")
            raise

    @rpc #(flags={"async"})
    def update_variable(self, variable_name, value):
        """
        update an ExperimentVariable by setting it to value

        this avoid updating the dataset, which would change the default on the next experiment we run
        and could therefore lead to unexpected behavior.
        """
        setattr(self, variable_name, value)

    @kernel
    def hardware_init(self):
        self.base.initialize_hardware()

    def run(self):

        for variable1_value in self.scan_sequence1:
            setattr(self, self.scan_variable1, variable1_value)
            self.hardware_init()

            for measurement in range(self.n_measurements):
                self.experiment_function()

    @kernel
    def initialize_datasets(self):
        self.set_dataset("photocounts", [0])
        self.set_dataset("photocounts2", [0])
        self.set_dataset("photocount_bins", [50], broadcast=True)

    # @kernel
    # def reinitialize_datasets(self):
    #     # todo
    #
    @kernel
    def update_datasets(self,variable1_value):
        self.set_dataset(self.scan_variable1, variable1_value) # don't broadcast

    @kernel
    def experiment_loop(self):
        """
        The experiment loop.

        Anything that could be experiment-specific, such as photocount datasets, should be dealt
        with in experiment_waveform.
        :return:
        """

        # todo: set datasets that will be useful for plotting

        # for variable1_value in self.scan_sequence1:

        # # this doesn't succeed at updating the variable.
        # self.update_variable(self.scan_variable1, variable1_value)
        # # self.set_dataset(self.scan_variable1, variable1_value) # don't broadcast
        #
        # # for measurement in ...
        self.experiment_function(self)


        self.print_async("Experiment loop finished")



