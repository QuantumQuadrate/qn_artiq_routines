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

# this is where your experiment function should live
from subroutines.experiment_functions import *
import subroutines.experiment_functions as exp_functions

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

        experiment_function_names_list = [x for x in dir(exp_functions)
            if ('__' not in x and str(type(getattr(exp_functions,x)))=="<class 'function'>"
                and 'experiment' in x)]

        # a function that take no arguments that gets imported and run
        # self.setattr_argument('experiment_function', StringValue('test'))
        self.setattr_argument('experiment_function',
                              EnumerationValue(experiment_function_names_list))

        # toggles an interleaved control experiment, but what this means or whether
        # it has an effect depends on experiment_function
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
        self.scan_sequence1 = eval(self.scan_sequence1)
        self.n_iterations1 = len(self.scan_sequence1) # this might not get used

        if self.scan_variable2 != '':
            assert hasattr(self, self.scan_variable2), (f"There is no ExperimentVariable " + self.scan_variable2 +
                                                        ". Did you mistype it?")
            self.scan_sequence2 = eval(self.scan_sequence2)
            # todo: add some error handling
            self.n_iterations2 = len(self.scan_sequence2)
        else:
            self.scan_variable2 = None
            self.scan_sequence2 = np.zeros(1) # we need a non-empty array to have a definite type for ARTIQ
            self.n_iterations2 = 0

        try:
            self.experiment_name = self.experiment_function
            self.experiment_function = lambda :eval(self.experiment_name)(self)
        except NameError as e:
            print(f"The function {experiment_name} is not defined. Did you forget to import it?")
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

    def reset_datasets(self):
        """
        set datasets that are redefined each iteration.

        typically these datasets are used for plotting which would be meaningless if we continued to append to the data,
        e.g. for the second readout histogram which we expect in general will change as experiment parameters induce
        different amount of atom loss.
        :return:
        """
        self.set_dataset("test_dataset", [0], broadcast=True)
        self.set_dataset('photocounts_current_iteration', [0], broadcast=True)
        self.set_dataset('photocounts2_current_iteration', [0], broadcast=True)

    def run(self):
        """
        Step through the variable values defined by the scan sequences and run the experiment function.

        Because the scan variables can be any ExperimentVariable, which includes values used to initialize
        hardware (e.g. a frequency for a dds channel), the hardware is reinitialized in each step of the
        variable scan, i.e., each iteration.
        """

        self.initialize_datasets()

        # scan in up to 2 dimensions. for each setting of the parameters, run experiment_function n_measurement times
        iteration = 0
        for variable1_value in self.scan_sequence1:
            # update the variable
            setattr(self, self.scan_variable1, variable1_value)

            for variable2_value in self.scan_sequence2:

                self.set_dataset("iteration", iteration, broadcast=True)

                if self.scan_variable2 != None:
                    setattr(self, self.scan_variable2, variable2_value)

                self.hardware_init()
                self.reset_datasets()

                # the measurement loop.
                self.experiment_function()

                iteration += 1

        self.write_results({'name':self.experiment_name[:-11]})  # write the h5 file here in case worker refuses to die





