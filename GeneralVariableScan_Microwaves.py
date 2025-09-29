"""
GeneralVariableScan for Microwave transition scans

compatible with GVS analysis

Added fitting functionality
"""


from artiq.experiment import *
from artiq.language import us, ns, MHz
import logging

import numpy as np
from numpy import array  # necessary for some override_ExperimentVariable entries

import sys, os

cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")
from utilities.BaseExperiment import BaseExperiment

# this is where your experiment function should live
from subroutines.experiment_functions import *
import subroutines.experiment_functions as exp_functions
from subroutines.aom_feedback import AOMPowerStabilizer


#### fitting functions
from fitting.rabi_flop import rabi_flop
from fitting.rabi_flop_reversed import rabi_flop_reversed
from fitting.resonance_dip import resonance_dip
from fitting.resonance_positive_dip import resonance_positive_dip

# from fitting_oxford.lorentzian import lorentzian
# from fitting.gaussian import gaussian


fit_model_dict = {
    "rabi_flop": rabi_flop,
    "rabi_flop_reversed": rabi_flop_reversed,
    "resonance_dip": resonance_dip,
    "resonance_positive_dip": resonance_positive_dip,
    # "gaussian": gaussian
}


class GeneralVariableScan_Microwaves(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        # the number of measurements to be made for a certain setting of the
        # experiment parameters
        self.setattr_argument("n_measurements", NumberValue(100, ndecimals=0, step=1), "General Setting")
        self.setattr_argument('override_ExperimentVariables', StringValue("{'dummy_variable':4}"), "General Setting")
        self.setattr_argument('enable_fitting', BooleanValue(default=False), "General Setting")
        self.setattr_argument('update_dataset', BooleanValue(default=False), "General Setting")

        self.setattr_argument("Frequency_00_Scan", BooleanValue(default=False),"Microwave Scans - choose one of the following")
        self.setattr_argument("Frequency_01_Scan", BooleanValue(default=False),"Microwave Scans - choose one of the following")
        self.setattr_argument("Frequency_11_Scan", BooleanValue(default=False),"Microwave Scans - choose one of the following")

        self.setattr_argument("Time_00_Scan", BooleanValue(default=False), "Microwave Scans - choose one of the following")
        self.setattr_argument("Time_01_Scan", BooleanValue(default=False), "Microwave Scans - choose one of the following")
        self.setattr_argument("Time_11_Scan", BooleanValue(default=False), "Microwave Scans - choose one of the following")

        self.setattr_argument("Ramsey_00_Scan", BooleanValue(default=False), "Microwave Scans - choose one of the following")
        self.setattr_argument("Ramsey_01_Scan", BooleanValue(default=False), "Microwave Scans - choose one of the following")
        self.setattr_argument("Ramsey_11_Scan", BooleanValue(default=False), "Microwave Scans - choose one of the following")

        #### Frequency scan variables
        self.setattr_argument("freq_scan_range_left_kHz", NumberValue(100.0, ndecimals=1, step=1), "Freq Scan variables - centered at resonance")
        self.setattr_argument("freq_scan_range_right_kHz", NumberValue(100.0, ndecimals=1, step=1), "Freq Scan variables - centered at resonance")
        self.setattr_argument("freq_scan_step_size_kHz", NumberValue(20.0, ndecimals=1, step=1), "Freq Scan variables - centered at resonance")

        #### Time scan variables
        self.setattr_argument("time_scan_range_start_us", NumberValue(0.0, ndecimals=2, step=1), "Time Scan variables")
        self.setattr_argument("time_scan_range_end_us", NumberValue(15.0, ndecimals=2, step=1), "Time Scan variables")
        self.setattr_argument("time_scan_step_size_us", NumberValue(1.0, ndecimals=2, step=1), "Time Scan variables")

        self.base.set_datasets_from_gui_args()
        print("build - done")


    def prepare(self):
        """
        performs initial calculations and sets parameter values before
        running the experiment.
        """
        self.base.prepare()

        self.override_ExperimentVariables_dict = eval(self.override_ExperimentVariables)
        assert type(self.override_ExperimentVariables_dict) == dict, \
            "override_ExperimentVariables should be a python dictionary"

        if self.Frequency_00_Scan:
            print("Frequency_00_Scan with pi pulse")
            self.override_ExperimentVariables_dict["t_microwave_pulse"] = self.t_microwave_00_pulse

            self.scan_variable1_name = 'f_microwaves_dds'
            self.scan_variable1 = str(self.scan_variable1_name)

            self.scan_sequence1 = (np.arange(self.f_microwaves_00_dds - self.freq_scan_range_left_kHz * kHz,
                                                 self.f_microwaves_00_dds + self.freq_scan_range_right_kHz * kHz,
                                       self.freq_scan_step_size_kHz * kHz))

            ### experiment function
            if self.which_node == 'bob':
                self.experiment_name = "microwave_Rabi_2_CW_OP_UW_FORT_experiment"
            elif self.which_node == 'alice':
                self.experiment_name = "microwave_Rabi_2_experiment"

            fit_model = "resonance_dip"
            self.which_fit_model = fit_model_dict[fit_model]

        elif self.Frequency_01_Scan:
            print("Frequency_01_Scan with pi pulse")
            self.override_ExperimentVariables_dict["t_microwave_pulse"] = self.t_microwave_01_pulse

            self.scan_variable1_name = 'f_microwaves_dds'
            self.scan_variable1 = str(self.scan_variable1_name)
            self.scan_sequence1 = (np.arange(self.f_microwaves_01_dds - self.freq_scan_range_left_kHz * kHz,
                                                 self.f_microwaves_01_dds + self.freq_scan_range_right_kHz * kHz,
                                       self.freq_scan_step_size_kHz * kHz))

            ### experiment function
            if self.which_node == 'bob':
                self.experiment_name = "microwave_Rabi_2_CW_OP_UW_FORT_experiment"
            elif self.which_node == 'alice':
                self.experiment_name = "microwave_Rabi_2_experiment"

            fit_model = "resonance_dip"
            self.which_fit_model = fit_model_dict[fit_model]

        elif self.Frequency_11_Scan:
            print("Frequency_11_Scan with pi pulse")
            self.override_ExperimentVariables_dict["t_microwave_11_pulse"] = self.t_microwave_11_pulse

            self.scan_variable1_name = 'f_microwaves_11_dds'
            self.scan_variable1 = str(self.scan_variable1_name)

            self.scan_sequence1 = (np.arange(self.f_microwaves_11_dds - self.freq_scan_range_left_kHz * kHz,
                                                 self.f_microwaves_11_dds + self.freq_scan_range_right_kHz * kHz,
                                       self.freq_scan_step_size_kHz * kHz))

            ### experiment function
            if self.which_node == 'bob':
                self.experiment_name = "microwave_Rabi_2_CW_OP_UW_FORT_11_experiment"
            elif self.which_node == 'alice':
                self.experiment_name = "microwave_Rabi_2_experiment"

            fit_model = "resonance_positive_dip"
            self.which_fit_model = fit_model_dict[fit_model]


        elif self.Time_00_Scan:
            print("Time_00_Scan")
            self.override_ExperimentVariables_dict["f_microwaves_dds"] = self.f_microwaves_00_dds

            self.scan_variable1_name = 't_microwave_pulse'
            self.scan_variable1 = str(self.scan_variable1_name)

            self.scan_sequence1 = (np.arange(self.time_scan_range_start_us * us,
                                             self.time_scan_range_end_us * us,
                                             self.time_scan_step_size_us * us))

            ### experiment function
            if self.which_node == 'bob':
                self.experiment_name = "microwave_Rabi_2_CW_OP_UW_FORT_experiment"
            elif self.which_node == 'alice':
                self.experiment_name = "microwave_Rabi_2_experiment"

            fit_model = "rabi_flop"
            self.which_fit_model = fit_model_dict[fit_model]

        elif self.Time_01_Scan:
            print("Time_01_Scan")
            self.override_ExperimentVariables_dict["f_microwaves_dds"] = self.f_microwaves_01_dds

            self.scan_variable1_name = 't_microwave_pulse'
            self.scan_variable1 = str(self.scan_variable1_name)

            self.scan_sequence1 = (np.arange(self.time_scan_range_start_us * us,
                                             self.time_scan_range_end_us * us,
                                             self.time_scan_step_size_us * us))

            ### experiment function
            if self.which_node == 'bob':
                self.experiment_name = "microwave_Rabi_2_CW_OP_UW_FORT_experiment"
            elif self.which_node == 'alice':
                self.experiment_name = "microwave_Rabi_2_experiment"

            fit_model = "rabi_flop"
            self.which_fit_model = fit_model_dict[fit_model]

        elif self.Time_11_Scan:
            print("Time_11_Scan")
            self.override_ExperimentVariables_dict["f_microwaves_11_dds"] = self.f_microwaves_11_dds

            self.scan_variable1_name = 't_microwave_11_pulse'
            self.scan_variable1 = str(self.scan_variable1_name)

            self.scan_sequence1 = (np.arange(self.time_scan_range_start_us * us,
                                             self.time_scan_range_end_us * us,
                                             self.time_scan_step_size_us * us))

            ### experiment function
            if self.which_node == 'bob':
                self.experiment_name = "microwave_Rabi_2_CW_OP_UW_FORT_11_experiment"
            elif self.which_node == 'alice':
                self.experiment_name = "microwave_Rabi_2_experiment"

            fit_model = "rabi_flop_reversed"
            # fit_model = "rabi_flop"
            self.which_fit_model = fit_model_dict[fit_model]

        elif self.Ramsey_00_Scan:
            print("Ramsey_00_Scan")
            self.override_ExperimentVariables_dict["f_microwaves_dds"] = self.f_microwaves_00_dds
            self.override_ExperimentVariables_dict["t_microwave_pulse"] = self.t_microwave_00_pulse

            self.scan_variable1_name = 't_delay_between_shots'
            self.scan_variable1 = str(self.scan_variable1_name)

            self.scan_sequence1 = (np.arange(self.time_scan_range_start_us * us,
                                             self.time_scan_range_end_us * us,
                                             self.time_scan_step_size_us * us))

            ### experiment function
            if self.which_node == 'bob':
                self.experiment_name = "microwave_Rabi_2_CW_OP_UW_FORT_Ramsey_experiment"
            elif self.which_node == 'alice':
                self.experiment_name = "microwave_Rabi_2_experiment"

            ###todo:  fitting not ready for Ramsey yet
            self.enable_fitting = False

        elif self.Ramsey_01_Scan:
            print("Ramsey_01_Scan")
            self.override_ExperimentVariables_dict["f_microwaves_dds"] = self.f_microwaves_01_dds
            self.override_ExperimentVariables_dict["t_microwave_pulse"] = self.t_microwave_01_pulse

            self.scan_variable1_name = 't_delay_between_shots'
            self.scan_variable1 = str(self.scan_variable1_name)

            self.scan_sequence1 = (np.arange(self.time_scan_range_start_us * us,
                                             self.time_scan_range_end_us * us,
                                             self.time_scan_step_size_us * us))

            ### experiment function
            if self.which_node == 'bob':
                self.experiment_name = "microwave_Rabi_2_CW_OP_UW_FORT_Ramsey_experiment"
            elif self.which_node == 'alice':
                self.experiment_name = "microwave_Rabi_2_experiment"

            ###todo:  fitting not ready for Ramsey yet
            self.enable_fitting = False

        elif self.Ramsey_11_Scan:
            print("Ramsey_11_Scan")

            self.scan_variable1_name = 't_delay_between_shots'
            self.scan_variable1 = str(self.scan_variable1_name)

            self.scan_sequence1 = (np.arange(self.time_scan_range_start_us * us,
                                             self.time_scan_range_end_us * us,
                                             self.time_scan_step_size_us * us))

            ### experiment function
            if self.which_node == 'bob':
                self.experiment_name = "microwave_Rabi_2_CW_OP_UW_FORT_11_Ramsey_experiment"
            elif self.which_node == 'alice':
                self.experiment_name = "microwave_Rabi_2_experiment"

            ###todo:  fitting not ready for Ramsey yet
            self.enable_fitting = False

        else:
            assert "Microwave Scan is not selected!"

        #### some more lines to be compatible with GVS analysis
        self.experiment_function = lambda: eval(self.experiment_name)(self)

        self.n_iterations1 = len(self.scan_sequence1) # this might not get used
        self.scan_variable2_name = ''
        self.scan_variable2 = None
        self.scan_sequence2 = np.zeros(1) # we need a non-empty array to have a definite type for ARTIQ
        self.n_iterations2 = 0

        scan_vars = [self.scan_variable1_name, self.scan_variable2_name]
        scan_vars = [x for x in scan_vars if x != '']
        self.scan_var_labels = ','.join(scan_vars)
        self.scan_var_filesuffix = '_and_'.join(scan_vars)

        self.set_dataset('scan_variable1_name', self.scan_variable1_name, broadcast=True)
        self.set_dataset('scan_sequence1', self.scan_sequence1, broadcast=True)
        self.set_dataset('scan_variable2_name', self.scan_variable2_name, broadcast=True)
        self.set_dataset('scan_sequence2', self.scan_sequence2, broadcast=True)

        self.set_dataset('experiment_function', self.experiment_name, broadcast=True)


        for variable, value in self.override_ExperimentVariables_dict.items():
            assert hasattr(self, variable), (f"There is no ExperimentVariable " + variable +
                                                    ". Did you mistype it?")

        self.measurement = 0
        self.SPCM0_RO1 = 0
        self.SPCM0_RO2 = 0
        self.SPCM1_RO1 = 0
        self.SPCM1_RO2 = 0

        # if there are multiple experiments in the schedule, then there might be something that has updated the datasets
        # e.g., as a result of an optimization scan. We want to make sure that this experiment uses the most up-to-date
        # datasets. However, ARTIQ runs build and prepare while the previous experiment is running, so our base.build
        # and base.prepare calls might use outdated datasets. This checks if there are other experiments scheduled
        # earlier than this one, and sets a flag so that base.build and base.prepare can be called again if needed.
        status_dict = self.scheduler.get_status()
        my_rid = self.scheduler.rid
        earlier_experiments = len([rid for rid, _ in status_dict.items() if rid < my_rid])
        logging.info("my rid is", my_rid, ", and there are", earlier_experiments, " experiment(s) that I am waiting on to run")
        self.needs_fresh_build = earlier_experiments > 0

    @kernel
    def initialize_hardware(self):
        self.base.initialize_hardware()

    def initialize_datasets(self):
        self.base.initialize_datasets()

        self.set_dataset(self.scan_var_dataset,self.scan_var_labels,broadcast=True)
        self.set_dataset(self.scan_sequence1_dataset,self.scan_sequence1, broadcast=True)
        self.set_dataset(self.scan_sequence2_dataset,self.scan_sequence2, broadcast=True)

        for var,val in self.override_ExperimentVariables_dict.items():
            self.set_dataset(var, val)

        value = 0.0
        for ch_i in range(len(self.laser_stabilizer.all_channels)):
            self.set_dataset(self.laser_stabilizer.all_channels[ch_i].dB_history_dataset,
                             [float(self.initial_RF_dB_values[ch_i])], broadcast=True, persist=True)

    def reset_datasets(self):
        """
        set datasets that are redefined each iteration.

        typically these datasets are used for plotting which would be meaningless if we continued to append to the data,
        e.g. for the second readout histogram which we expect in general will change as experiment parameters induce
        different amount of atom loss.
        :return:
        """
        self.set_dataset("test_dataset", [0], broadcast=True)
        self.set_dataset('SPCM0_RO1_current_iteration', [0], broadcast=True)
        self.set_dataset('SPCM0_RO2_current_iteration', [0], broadcast=True)
        self.set_dataset('SPCM1_RO1_current_iteration', [0], broadcast=True)
        self.set_dataset('SPCM1_RO2_current_iteration', [0], broadcast=True)
        self.set_dataset('BothSPCMs_RO1_current_iteration', [0], broadcast=True)
        self.set_dataset('BothSPCMs_RO2_current_iteration', [0], broadcast=True)

        # these are set here because running BaseExperiment.initialize_hardware resets these to be empty
        self.set_dataset(self.scan_var_dataset, self.scan_var_labels, broadcast=True)
        self.set_dataset(self.scan_sequence1_dataset, self.scan_sequence1, broadcast=True)
        self.set_dataset(self.scan_sequence2_dataset, self.scan_sequence2, broadcast=True)

    def initialize_dependent_variables(self):
        """
        anything that happens in base.prepare, but can add other things as well

        this is what allows us to adjust the setpoints for the AOMStabilizer, which needs to be
        reinstantiated before we call the experiment function.
        """
        self.base.prepare()

    @kernel
    def warm_up(self):
        """hardware init and turn things on"""

        self.core.reset()

        self.zotino0.set_dac(
            [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
            channels=self.coil_channels)

        self.dds_cooling_DP.sw.on()
        self.dds_AOM_A1.sw.on()
        self.dds_AOM_A2.sw.on()
        self.dds_AOM_A3.sw.on()
        self.dds_AOM_A4.sw.on()
        self.dds_AOM_A5.sw.on()
        self.dds_AOM_A6.sw.on()
        self.dds_FORT.sw.on()

        # delay for AOMs to thermalize
        delay(200 * ms)

        self.core.break_realtime()
        # warm up to get make sure we get to the setpoints

        if self.enable_laser_feedback:
            for i in range(10):
                self.laser_stabilizer.run()
        else:
            delay(200*ms) # lotsa slack
        self.dds_FORT.sw.on()

    def run(self):
        """
        Step through the variable values defined by the scan sequences and run the experiment function.

        Because the scan variables can be any ExperimentVariable, which includes values used to initialize
        hardware (e.g. a frequency for a dds channel), the hardware is reinitialized in each step of the
        variable scan, i.e., each iteration.
        """

        if self.needs_fresh_build:
            self.base.build()
            self.base.prepare()

        self.initialize_datasets()

        # scan in up to 2 dimensions. for each setting of the parameters, run experiment_function n_measurement times
        iteration = 0
        self.set_dataset("iteration", iteration, broadcast=True)

        # override specific variables. this will apply to the entire scan, so it is outside the loops
        for variable, value in self.override_ExperimentVariables_dict.items():
            setattr(self, variable, value)

        # self.warm_up()

        ##### Scan sequence
        for variable1_value in self.scan_sequence1:
            # update the variable. setattr can't be called on the kernel, and this is what
            # allows us to update an experiment variable without hardcoding it, i.e.
            # explicitly naming the variable. that is why this run method does not
            # have a kernel decorator, and we have to re-initialize the hardware each
            # iteration.
            setattr(self, self.scan_variable1, variable1_value)
            logging.info(f"current iteration: {self.scan_variable1_name} = {variable1_value}")

            for variable2_value in self.scan_sequence2:

                self.set_dataset("iteration", iteration, broadcast=True)

                if self.scan_variable2 != None:
                    setattr(self, self.scan_variable2, variable2_value)
                    logging.info(f"current iteration: {self.scan_variable2_name} ={variable2_value}")

                self.initialize_dependent_variables()
                self.initialize_hardware()
                self.reset_datasets()

                # the measurement loop.
                self.experiment_function()

                # write and overwrite the file here so we can quit the experiment early without losing data
                self.write_results({'name': self.experiment_name[:-11] + "_scan_over_" + self.scan_var_filesuffix})

                iteration += 1

        #### for fitting
        if self.enable_fitting:
            BothSPCMs_RO1 = self.get_dataset("BothSPCMs_RO1")
            BothSPCMs_RO2 = self.get_dataset("BothSPCMs_RO2")

            retention_array = self.get_retention(BothSPCMs_RO1, BothSPCMs_RO2, self.n_measurements, len(self.scan_sequence1),
                               self.single_atom_threshold * self.t_SPCM_first_shot)



            #  ### todo: this is for resonance_dip: 00 transition:  f0 = 334.662MHz
            # self.scan_sequence1 = np.array([3.34565e+08, 3.34575e+08, 3.34585e+08, 3.34595e+08, 3.34605e+08,
            # 3.34615e+08, 3.34625e+08, 3.34635e+08, 3.34645e+08, 3.34655e+08,
            # 3.34665e+08, 3.34675e+08, 3.34685e+08, 3.34695e+08, 3.34705e+08,
            # 3.34715e+08, 3.34725e+08, 3.34735e+08, 3.34745e+08, 3.34755e+08])
            #
            # retention_array = np.array([0.4893617 , 0.42777778, 0.34078212, 0.3038674 , 0.27071823,
            # 0.15384615, 0.16292135, 0.05586592, 0.08333333, 0.08839779,
            # 0.04945055, 0.06741573, 0.10169492, 0.11666667, 0.17127072,
            # 0.21590909, 0.27513228, 0.34920635, 0.37837838, 0.43243243])
            #
            # ### todo: this is for resonance_dip: 01 transition:  f0 = 337.322MHz
            # self.scan_sequence1 = np.array([3.37209265e+08, 3.37219265e+08, 3.37229265e+08, 3.37239265e+08,
            #        3.37249265e+08, 3.37259265e+08, 3.37269265e+08, 3.37279265e+08,
            #        3.37289265e+08, 3.37299265e+08, 3.37309265e+08, 3.37319265e+08,
            #        3.37329265e+08, 3.37339265e+08, 3.37349265e+08, 3.37359265e+08,
            #        3.37369265e+08, 3.37379265e+08, 3.37389265e+08, 3.37399265e+08])
            #
            # retention_array = np.array([0.96354167, 0.9025641 , 0.84615385, 0.80526316, 0.65803109,
            #        0.55026455, 0.4507772 , 0.34408602, 0.21195652, 0.13089005,
            #        0.08152174, 0.03763441, 0.04787234, 0.08695652, 0.1657754 ,
            #        0.24736842, 0.35978836, 0.50810811, 0.6       , 0.75520833])

            # ### todo: fit this to resonance_positive_dip: 11 transition: f0 = 339.993
            # self.scan_sequence1 = np.array([3.39843e+08, 3.39853e+08, 3.39863e+08, 3.39873e+08, 3.39883e+08,
            #         3.39893e+08, 3.39903e+08, 3.39913e+08, 3.39923e+08, 3.39933e+08,
            #         3.39943e+08, 3.39953e+08, 3.39963e+08, 3.39973e+08, 3.39983e+08,
            #         3.39993e+08, 3.40003e+08, 3.40013e+08, 3.40023e+08, 3.40033e+08])
            #
            # retention_array = np.array([0.23783784, 0.26666667, 0.32022472, 0.40935673, 0.43455497,
            #         0.52150538, 0.59782609, 0.60526316, 0.68108108, 0.72727273,
            #         0.8       , 0.83246073, 0.9144385 , 0.93406593, 0.96315789,
            #         0.95721925, 0.91666667, 0.92553191, 0.88947368, 0.86263736])

            # ### todo: time scan for 00 transition (4.05us)
            # self.scan_sequence1 = np.array([0.e+00, 1.e-06, 2.e-06, 3.e-06, 4.e-06, 5.e-06, 6.e-06, 7.e-06])
            #
            # retention_array = np.array([0.95876289, 0.82446809, 0.51308901, 0.19047619, 0.08695652,
            #        0.15, 0.44565217, 0.74093264])

            #  ### todo: fit this to reversed rabi (11 transition )- t_pi = 4.4us
            #      self.scan_sequence1 = np.array([0.e+00, 1.e-06, 2.e-06, 3.e-06, 4.e-06, 5.e-06, 6.e-06, 7.e-06,
            # 8.e-06, 9.e-06])
            #      retention_array = np.array([0.16489362, 0.18918919, 0.51578947, 0.69948187, 0.9119171 ,
            # 0.82010582, 0.70491803, 0.42408377, 0.19889503, 0.16201117])

            t = self.scan_sequence1
            y = retention_array

            p, p_err = self.which_fit_model.fit(t, y, evaluate_function=False)
            print(p)

            # print(p["f0"])
            # print(type(p["f0"]))

            if self.update_dataset:
                if self.Frequency_00_Scan:
                    self.set_dataset("f_microwaves_00_dds", round(float(p["f0"]), 0), broadcast=True, persist=True)
                    print(self.scan_variable1_name, " updated to ", p["f0"])

                elif self.Frequency_01_Scan:
                    self.set_dataset("f_microwaves_01_dds", round(float(p["f0"]), 0), broadcast=True, persist=True)
                    print(self.scan_variable1_name, " updated to ", round(float(p["f0"]), 1))

                elif self.Frequency_11_Scan:
                    if float(p["sigma"]) > 130e3:
                        print("fitted sigma > 130kHz - 01 transition frequency retune seems to be necessary")
                    self.set_dataset("f_microwaves_11_dds", round(float(p["f0"]), 0), broadcast=True, persist=True)
                    print(self.scan_variable1_name, " updated to ", p["f0"])

                elif self.Time_00_Scan:
                    self.set_dataset("t_microwave_00_pulse", round(float(p["t_pi"]),8), broadcast=True, persist=True)
                    print("t_microwave_00_pulse updated to ", round(float(p["t_pi"]),8))

                elif self.Time_01_Scan:
                    self.set_dataset("t_microwave_01_pulse", round(float(p["t_pi"]),8), broadcast=True, persist=True)
                    print("t_microwave_01_pulse updated to ", round(float(p["t_pi"]),8))

                elif self.Time_11_Scan:
                    self.set_dataset("t_microwave_11_pulse", round(float(p["t_pi"]), 8), broadcast=True, persist=True)
                    print("t_microwave_11_pulse updated to ", round(float(p["t_pi"]), 8))

            # round(float(p["f0"]), 1)



        print("****************    General Variable Scan DONE   *****************")


    def get_loading_and_retention(self, photocounts, photocounts2, measurements, iterations, cutoff):
        """
        Returns retention, loading rate, and number of atoms loaded for each experiment iteration.
        cutoff1 and cutoff2 (optional) are the atom loading thresholds in units counts.
        """

        retention_array = np.zeros(iterations)
        loading_rate_array = np.zeros(iterations)
        n_atoms_loaded_array = np.zeros(iterations)

        for i in range(iterations):
            shot1 = photocounts[i * measurements:(i + 1) * measurements]
            shot2 = photocounts2[i * measurements:(i + 1) * measurements]

            atoms_loaded = [x > cutoff for x in shot1]
            n_atoms_loaded = sum(atoms_loaded)
            atoms_retained = [x > cutoff and y for x, y in zip(shot2, atoms_loaded)]
            retention_fraction = 0 if not n_atoms_loaded > 0 else sum(atoms_retained) / sum(atoms_loaded)
            loading_rate_array[i] = n_atoms_loaded / measurements
            n_atoms_loaded_array[i] = n_atoms_loaded
            retention_array[i] = retention_fraction
        return retention_array, loading_rate_array, n_atoms_loaded_array

    def get_retention(self, photocounts, photocounts2, measurements, iterations, cutoff):
        """
        Returns retention, loading rate, and number of atoms loaded for each experiment iteration.
        cutoff1 and cutoff2 (optional) are the atom loading thresholds in units counts.
        """

        retention_array = np.zeros(iterations)

        for i in range(iterations):
            shot1 = photocounts[i * measurements:(i + 1) * measurements]
            shot2 = photocounts2[i * measurements:(i + 1) * measurements]

            atoms_loaded = [x > cutoff for x in shot1]
            n_atoms_loaded = sum(atoms_loaded)
            atoms_retained = [x > cutoff and y for x, y in zip(shot2, atoms_loaded)]
            retention_fraction = 0 if not n_atoms_loaded > 0 else sum(atoms_retained) / sum(atoms_loaded)
            retention_array[i] = retention_fraction

        return retention_array

