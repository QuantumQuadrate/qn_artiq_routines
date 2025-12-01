"""
GeneralVariableScan for Microwave transition scans

- compatible with GVS analysis
- enable_fitting: enables fitting of the scanned results.
- update_dataset: updates the dataset to the optimized fitted value.
"""


from artiq.experiment import *
from artiq.language import us, ns, MHz
import logging

import numpy as np
from numpy import array  # necessary for some override_ExperimentVariable entries

import copy
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
from fitting.resonance_peak import resonance_peak

# from fitting_oxford.lorentzian import lorentzian
# from fitting.gaussian import gaussian


fit_model_dict = {
    "rabi_flop": rabi_flop,
    "rabi_flop_reversed": rabi_flop_reversed,
    "resonance_dip": resonance_dip,
    "resonance_peak": resonance_peak,
    # "gaussian": gaussian
}

scan_options = [
    "Frequency_00_Scan", "Frequency_01_Scan", "Frequency_11_Scan",
    "Frequency_m10_Scan", "Frequency_m11_Scan",
    "Time_00_Scan", "Time_01_Scan", "Time_11_Scan", "Time_m10_Scan",
    "Ramsey_00_Scan", "Ramsey_01_Scan", "Ramsey_11_Scan",
]

scan_dict={
    "Frequency_00_Scan":{
        "print_statement": "Frequency_00_Scan with pi pulse",
        "override_items": {
            "t_microwave_pulse": "t_microwave_00_pulse"
        },

        "scan_variable1_name": "f_microwaves_dds",
        "center": "f_microwaves_00_dds",

        "experiment_name_alice": "microwave_Rabi_2_experiment",
        "experiment_name_bob": "microwave_Rabi_2_CW_OP_UW_FORT_experiment",

        "fit_model": "resonance_dip",
        "initialise": {
            "sigma": 80e3
        },
        "health_check_dataset_name": "health_check_uw_freq00"
    },

    "Frequency_01_Scan": {
        "print_statement": "Frequency_01_Scan with pi pulse",
        "override_items": {
            "t_microwave_pulse": "t_microwave_01_pulse"
        },

        "scan_variable1_name": "f_microwaves_dds",
        "center": "f_microwaves_01_dds",

        "experiment_name_alice": "microwave_Rabi_2_experiment",
        "experiment_name_bob": "microwave_Rabi_2_CW_OP_UW_FORT_experiment",

        "fit_model": "resonance_dip",
        "initialise": {
            "sigma": 50e3
        },
        "health_check_dataset_name": "health_check_uw_freq01"
    },

    "Frequency_11_Scan": {
        "print_statement": "Frequency_11_Scan with pi pulse",
        "override_items": {},

        "scan_variable1_name": "f_microwaves_11_dds",
        "center": "f_microwaves_11_dds",

        "experiment_name_alice": "microwave_map01_map11_experiment",
        "experiment_name_bob": "microwave_Rabi_2_CW_OP_UW_FORT_11_experiment",

        "fit_model": "resonance_peak",
        "initialise": {
            "sigma": 80e3
        },
        "health_check_dataset_name": "health_check_uw_freq11"
    },

    "Frequency_m10_Scan": {
        "print_statement": "Frequency_10_Scan with pi pulse",
        "override_items": {},

        "scan_variable1_name": "f_microwaves_m10_dds",
        "center": "f_microwaves_m10_dds",

        "experiment_name_alice": "microwave_map00_map0m1_experiment",
        "experiment_name_bob": "microwave_Rabi_2_CW_OP_UW_FORT_m10_experiment",

        "fit_model": "resonance_peak",
        "initialise": {},
        "health_check_dataset_name": "health_check_uw_freqm10"
    },

    "Frequency_m11_Scan": {
        "print_statement": "Frequency_m11_Scan with pi pulse",
        "override_items": {
            "t_MW_RF_pulse": "t_MW_RF_pulse"
        },

        "scan_variable1_name": "f_microwaves_m11_dds",
        "center": "f_microwaves_m11_dds",

        "experiment_name_alice": "microwave_map01_MWRFm11_experiment",
        "experiment_name_bob": "microwave_Rabi_2_CW_OP_UW_FORT_m10_experiment",

        "fit_model": "resonance_peak",
        "initialise": {},
        "health_check_dataset_name": "health_check_uw_freqm11"
    },

    # Time scans do not have health check dataset
    "Time_00_Scan": {
        "print_statement": "Time_00_Scan with freq00",
        "override_items": {
            "f_microwaves_dds": "f_microwaves_00_dds"
        },

        "scan_variable1_name": "t_microwave_pulse",
        "pi_pulse": "t_microwave_00_pulse",

        "experiment_name_alice": "microwave_Rabi_2_experiment",
        "experiment_name_bob": "microwave_Rabi_2_CW_OP_UW_FORT_experiment",

        "fit_model": "rabi_flop",
        "initialise": {},
        "health_check_dataset_name": "health_check_uw_freq00"
    },

    # Time scans do not have health check dataset
    "Time_01_Scan": {
        "print_statement": "Time_01_Scan with freq01",
        "override_items": {
            "f_microwaves_dds": "f_microwaves_01_dds"
        },

        "scan_variable1_name": "t_microwave_pulse",
        "pi_pulse": "t_microwave_01_pulse",

        "experiment_name_alice": "microwave_Rabi_2_experiment",
        "experiment_name_bob": "microwave_Rabi_2_CW_OP_UW_FORT_experiment",

        "fit_model": "rabi_flop",
        "initialise": {},
        "health_check_dataset_name": "health_check_uw_freq01"
    },

    # Time scans do not have health check dataset
    "Time_11_Scan": {
        "print_statement": "Time_11_Scan with freq11",
        "override_items": {},

        "scan_variable1_name": "t_microwave_11_pulse",
        "pi_pulse": "t_microwave_11_pulse",

        "experiment_name_alice": "microwave_map01_map11_experiment",
        "experiment_name_bob": "microwave_Rabi_2_CW_OP_UW_FORT_11_experiment",

        "fit_model": "rabi_flop_reversed",
        "initialise": {},
        "health_check_dataset_name": "health_check_uw_freq11"
    },

    # Time scans do not have health check dataset
    "Time_m10_Scan": {
        "print_statement": "Time_m10_Scan with freq m10",
        "override_items": {},

        "scan_variable1_name": "t_microwave_m10_pulse",
        "pi_pulse": "t_microwave_m10_pulse",

        "experiment_name_alice": "microwave_map00_map0m1_experiment",
        "experiment_name_bob": "microwave_Rabi_2_CW_OP_UW_FORT_m10_experiment",

        "fit_model": "rabi_flop_reversed",
        "initialise": {},

        "health_check_dataset_name": "health_check_uw_freqm10"
    },

    # Time scans do not have health check dataset
    "Ramsey_00_Scan": {
        "print_statement": "Ramsey_00_Scan",
        "override_items": {
            "f_microwaves_dds": "f_microwaves_00_dds",
            "t_microwave_pulse": "t_microwave_00_pulse"
        },

        "scan_variable1_name": "t_delay_between_shots",

        "experiment_name_alice": "microwave_Rabi_2_experiment",
        "experiment_name_bob": "microwave_Rabi_2_CW_OP_UW_FORT_Ramsey_experiment",

        ###todo:  fitting not ready for Ramsey yet
        "fit_model": "rabi_flop",
        "initialise": {}
    },

    # Time scans do not have health check dataset
    "Ramsey_01_Scan": {
        "print_statement": "Ramsey_01_Scan",
        "override_items": {
            "f_microwaves_dds": "f_microwaves_01_dds",
            "t_microwave_pulse": "t_microwave_01_pulse"
        },

        "scan_variable1_name": "t_delay_between_shots",

        "experiment_name_alice": "microwave_Rabi_2_experiment",
        "experiment_name_bob": "microwave_Rabi_2_CW_OP_UW_FORT_Ramsey_experiment",

        ###todo:  fitting not ready for Ramsey yet
        "fit_model": "rabi_flop",
        "initialise": {}
    },

    # Time scans do not have health check dataset
    "Ramsey_11_Scan": {
        "print_statement": "Ramsey_11_Scan",
        "override_items": {
        },

        "scan_variable1_name": "t_delay_between_shots",

        "experiment_name_alice": "microwave_Rabi_2_experiment",
        "experiment_name_bob": "microwave_Rabi_2_CW_OP_UW_FORT_11_Ramsey_experiment",

        ###todo:  fitting not ready for Ramsey yet
        "fit_model": "rabi_flop",
        "initialise": {}
    },

}




class MicrowaveScanOptimizer(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        # the number of measurements to be made for a certain setting of the
        # experiment parameters
        #todo: change this name to "run_health_check_and_schedule"
        #todo: add "schedule_multiple_scans"
        self.setattr_argument('run_health_check_and_optimize', BooleanValue(default=True), "Health Check")
        self.setattr_argument("target_fidelity", NumberValue(0.80, ndecimals=2, step=1), "Health Check")
        # this can be retrieved in  hdf5 file.

        self.setattr_argument("n_measurements", NumberValue(100, ndecimals=0, step=1), "General Scan Setting")
        self.setattr_argument('override_ExperimentVariables', StringValue("{'dummy_variable':4}"), "General Scan Setting")
        self.setattr_argument('enable_fitting', BooleanValue(default=True), "General Scan Setting")
        self.setattr_argument('enable_geometric_frequency_scan', BooleanValue(default=True), "General Scan Setting")

        self.setattr_argument("Frequency_00_Scan", BooleanValue(default=False),"Microwave Scans - Select one of the following")
        self.setattr_argument("Frequency_01_Scan", BooleanValue(default=False),"Microwave Scans - Select one of the following")
        self.setattr_argument("Frequency_11_Scan", BooleanValue(default=False),"Microwave Scans - Select one of the following")
        self.setattr_argument("Frequency_m10_Scan", BooleanValue(default=False),"Microwave Scans - Select one of the following")

        self.setattr_argument("Time_00_Scan", BooleanValue(default=False), "Microwave Scans - Select one of the following")
        self.setattr_argument("Time_01_Scan", BooleanValue(default=False), "Microwave Scans - Select one of the following")
        self.setattr_argument("Time_11_Scan", BooleanValue(default=False), "Microwave Scans - Select one of the following")
        self.setattr_argument("Time_m10_Scan", BooleanValue(default=False), "Microwave Scans - Select one of the following")

        self.setattr_argument("Ramsey_00_Scan", BooleanValue(default=False), "Microwave Scans - Select one of the following")
        self.setattr_argument("Ramsey_01_Scan", BooleanValue(default=False), "Microwave Scans - Select one of the following")
        self.setattr_argument("Ramsey_11_Scan", BooleanValue(default=False), "Microwave Scans - Select one of the following")

        #### Frequency scan variables
        self.setattr_argument("freq_scan_range_left_kHz", NumberValue(150.0, ndecimals=1, step=1), "[Default] Frequency Scan variables - centered at resonance")
        self.setattr_argument("freq_scan_range_right_kHz", NumberValue(150.0, ndecimals=1, step=1), "[Default] Frequency Scan variables - centered at resonance")
        self.setattr_argument("freq_scan_step_size_kHz", NumberValue(20.0, ndecimals=1, step=1), "[Default] Frequency Scan variables - centered at resonance")

        #### Frequency scan variables - for faster scan!
        self.setattr_argument("shrink_factor", NumberValue(2.5, ndecimals=1, step=1), "[Geometric] Frequency Scan variables - centered at resonance")
        self.setattr_argument("freq_scan_half_range_kHz", NumberValue(100.0, ndecimals=1, step=1), "[Geometric] Frequency Scan variables - centered at resonance")
        self.setattr_argument("freq_scan_min_step_size_kHz", NumberValue(10.0, ndecimals=1, step=1), "[Geometric] Frequency Scan variables - centered at resonance")

        #### Time scan variables
        self.setattr_argument("time_scan_sequence", StringValue('np.arange(0,10,1)*us'), "Time Scan variables")

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

        for variable, value in self.override_ExperimentVariables_dict.items():
            assert hasattr(self, variable), (f"There is no ExperimentVariable " + variable +
                                                    ". Did you mistype it?")

        # print(self.override_ExperimentVariables_dict)

        ### goes through the booleans and sets the scan type
        self.scan_type = self.get_scan_type()

        ### override items:
        for var, val in scan_dict[self.scan_type]["override_items"].items():
            self.override_ExperimentVariables_dict[var] = getattr(self, val)

        for variable, value in self.override_ExperimentVariables_dict.items():
            assert hasattr(self, variable), (f"There is no ExperimentVariable " + variable +
                                                    ". Did you mistype it?")

        # override specific variables. this will apply to the entire scan, so it is outside the loops
        for variable, value in self.override_ExperimentVariables_dict.items():
            setattr(self, variable, value)

        # print(self.override_ExperimentVariables_dict)

        ### setting the scan variable
        self.scan_variable1_name = scan_dict[self.scan_type]["scan_variable1_name"]
        self.scan_variable1 = str(self.scan_variable1_name)
        ### Note: there are some unnecessary lines so that the results are compatible with GVS analysis code

        ### setting the experiment function
        if self.which_node == 'bob':
            self.experiment_name = scan_dict[self.scan_type]["experiment_name_bob"]
        elif self.which_node == 'alice':
            self.experiment_name = scan_dict[self.scan_type]["experiment_name_alice"]
        self.experiment_function = lambda: eval(self.experiment_name)(self)

        ### setting up the fit model
        fit_model = scan_dict[self.scan_type]["fit_model"]
        self.which_fit_model = fit_model_dict[fit_model]
        self.initialise = scan_dict[self.scan_type]["initialise"]

        ### setting up the scan sequence
        if self.scan_type.startswith('Freq'):
            center = getattr(self, scan_dict[self.scan_type]["center"])
            print("center: ", center)
            if self.enable_geometric_frequency_scan:
                # todo: scan list must use mode = "pair", for compatibility with the adaptive scan logic
                self.scan_sequence1 = self.make_scan_list(center=center,
                                                          half_range_kHz=self.freq_scan_half_range_kHz,
                                                          min_step_kHz=self.freq_scan_min_step_size_kHz,
                                                          mode = "pair")
            else:
                self.scan_sequence1 = np.arange(center - self.freq_scan_range_left_kHz * kHz,
                                                center + self.freq_scan_range_right_kHz * kHz,
                                                self.freq_scan_step_size_kHz * kHz)

        elif self.scan_type.startswith(('Time', 'Ramsey')):
            self.scan_sequence1 = eval(self.time_scan_sequence)
        else:
            raise KeyError("There is no such thing as ", scan_variable1_name)


        ### Note: there are some unnecessary lines so that the results are compatible with GVS analysis code
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

        #todo: this only includes the experiments scheduled before
        #   ex1) When exp scheduled with due_date:
        #       if I scheduled an experiment and set the due date as tomorrow 7am.
        #       Then, I scheduled another experiment that starts now,
        #       the experiment that starts tomorrow "needs_fresh_build" but when it was scheduled,
        #       it did not have any earlier_experiments, so it will not build fresh.
        #   ex1 - sol) This can be solved by fresh building everytime when "due_date" != None
        #   ex1 - sol) When scheduled with due date, it starts to build when it is due.
        #              Thus, this problem does not occur.


        logging.info("my rid is", my_rid, ", and there are", earlier_experiments, " experiment(s) that I am waiting on to run")
        self.needs_fresh_build = earlier_experiments > 0

        # print(status_dict)

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
        print("self.n_measurements in GVS_Microwaves: ", self.n_measurements)
        print("self.shrink_factor in GVS_Microwaves: ", self.shrink_factor)

        # # override specific variables. this will apply to the entire scan, so it is outside the loops
        # for variable, value in self.override_ExperimentVariables_dict.items():
        #     setattr(self, variable, value)

        # todo: later just do all health checks no matter what
        if self.scan_type.startswith("Freq") and self.run_health_check_and_optimize:
            #todo: if there are more than 1 health checks, make a loop
            if self.health_check_general() == False:
                print("Initial Health Check - failed with fidelity: ", self.get_dataset(scan_dict[self.scan_type]["health_check_dataset_name"]))
                print("Scheduling Experiment for Optimization...")
                self.submit_opt_exp_general()
                # self.submit_opt_exp_general(override_arguments = {"freq_scan_half_range_kHz":150.0})

                self.scheduler.request_termination(self.scheduler.rid)

            self.core.comm.close()  # placing the hardware in a safe state and disconnecting it from the core device
            self.scheduler.pause()

            #todo: if running more than two health checks at0 the same time, fresh_build is necessary?
            #       => No. the things that are usually done in build/prepare is embedded in health_check
            #          Health Check only runs a single point scan. Thus, does not need more that that.


        else:
            ##### Scan sequence - same as GVS. (to be compatibel with original GVS analysis)
            ###todo: adaptive scan - make this visible in gui
            scan_with_fixed_sequence = False    ### adaptive scan enabled if False
            if scan_with_fixed_sequence:
                print("scanning with fixed sequence")
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

                        self.base.prepare()
                        self.initialize_hardware()
                        self.reset_datasets()

                        # the measurement loop.
                        self.experiment_function()

                        # write and overwrite the file here so we can quit the experiment early without losing data
                        self.write_results({'name': self.experiment_name[:-11] + "_scan_over_" + self.scan_var_filesuffix})

                        iteration += 1

            else:
                print("scanning with adpative sequence")
                #todo: scan list must use mode = "pair", for compatibility with the adaptive scan logic
                pending_sequence1 = list(self.scan_sequence1)
                did_refine = False

                scanned_points = []  # ← store every scanned variable1 value

                while pending_sequence1:
                    variable1_value = pending_sequence1.pop(0)

                    # ----- inner scanning loop -----
                    setattr(self, self.scan_variable1, variable1_value)

                    logging.info(f"current iteration: {self.scan_variable1_name} = {variable1_value}")

                    self.set_dataset("iteration", iteration, broadcast=True)

                    self.base.prepare()
                    self.initialize_hardware()
                    self.reset_datasets()

                    self.experiment_function()
                    self.write_results({
                        'name': self.experiment_name[:-11] + "_scan_over_" + self.scan_var_filesuffix
                    })

                    iteration += 1

                    # store the scanned variable1_value
                    scanned_points.append(variable1_value)

                    self.scan_sequence1 = np.array(scanned_points)
                    self.set_dataset("scan_sequence1", (scanned_points), broadcast=True)

                    if (not did_refine) and (len(scanned_points) == 4):
                        print("in decision loop")

                        BothSPCMs_RO1 = self.get_dataset("BothSPCMs_RO1")
                        BothSPCMs_RO2 = self.get_dataset("BothSPCMs_RO2")

                        retention_array, loading_rate_array, n_atoms_loaded_array = self.get_loading_and_retention(
                            BothSPCMs_RO1, BothSPCMs_RO2,
                            self.n_measurements,
                            int((len(BothSPCMs_RO1) - 1) / self.n_measurements),
                            self.single_atom_threshold * self.t_SPCM_first_shot
                        )

                        fit_model = scan_dict[self.scan_type]["fit_model"]
                        ### setting target retention to calculate fidelity
                        if fit_model in ["resonance_dip", "rabi_flop"]:
                            pass
                        elif fit_model in ["resonance_peak", "rabi_flop_reversed"]:
                            retention_array = 1.0 - retention_array


                        # First four points of the INITIAL scan sequence:
                        # index 0 -> -x, index 1 -> +x, index 2 -> -x/2, index 3 -> +x/2
                        f0, f1, f2, f3 = self.scan_sequence1[0:4]
                        R0, R1, R2, R3 = map(float, retention_array[0:4])

                        # slopes:
                        # left two points: (-x -> -x/2) => indices 0 -> 2
                        slope_left = R2 - R0
                        # right two points: (+x/2 -> +x) => indices 3 -> 1
                        slope_right = R1 - R3

                        lowest_index = int(np.argmin(retention_array[:4]))
                        lowest_value = float(retention_array[lowest_index])

                        threshold_low = 0.7 # tune this as needed
                        new_center = None

                        # 2.0 well centered:
                        if (lowest_index == 2 or lowest_index == 3) and slope_left < 0 and slope_right > 0:
                            print("Case 2.0: Well centered")
                            ### do something
                        # 2.1 left skewed:
                        # lowest at -x (index 0) and left slope positive
                        elif lowest_index == 0 and slope_left > 0 and lowest_value < threshold_low:
                            print("Case 2.1: left-skewed (min at -x)")
                            # example: shift center left by half current span
                            if lowest_value < 0.2:
                                new_center = f0 - self.freq_scan_min_step_size_kHz
                            elif lowest_value < 0.4:
                                step = abs(f2 - f0)  # distance between -x and -x/2
                                new_center = f0 - step  # move further left
                            else:
                                step = 2 * abs(f2 - f0)  # twice distance between -x and -x/2
                                new_center = f0 - step  # move further left

                        # 2.2 right skewed:
                        # lowest at +x (index 1) and right slope negative
                        elif lowest_index == 1 and slope_right < 0 and lowest_value < threshold_low:
                            print("Case 2.2: right-skewed (min at +x)")
                            if lowest_value < 0.2:
                                new_center = f1 + self.freq_scan_min_step_size_kHz
                            elif lowest_value < 0.4:
                                step = abs(f1 - f3)  # distance between +x/2 and +x
                                new_center = f1 + step  # move further right
                            else:
                                step = 2 * abs(f2 - f0)  # twice distance between -x and -x/2
                                new_center = f0 - step  # move further left

                        ##todo: case2.3 and 2.4 cases are mostly captured in 2.0 case. However, sometimes,
                        ## due to noise/uncertainty or due to sideband peak or etc slope condiction might not
                        ## satisfy case2.0; these lines are to prevent that
                        # 2.3 within range but skewed left:
                        # lowest at -x/2 (index 2) and slope_left negative
                        elif lowest_index == 2 and slope_left < 0 and lowest_value < threshold_low:
                            print("Case 2.3: in range but skewed left")
                            # e.g. pick a center between -x and +x/2 leaning left
                            if lowest_value < 0.2:
                                new_center = f2 + self.freq_scan_min_step_size_kHz
                            else:
                                step = abs((f2 - f0)/2)
                                new_center = f2 + step

                        # 2.4 within range but skewed right:
                        # lowest at +x/2 (index 3) and slope_right positive
                        elif lowest_index == 3 and slope_right > 0 and lowest_value < threshold_low:
                            print("Case 2.4: in range but skewed right")
                            if lowest_value < 0.2:
                                new_center = f3 - self.freq_scan_min_step_size_kHz
                            else:
                                step = abs((f1 - f3)/2)
                                new_center = f3 - step

                        else:
                            print("Not in case 2.1 ~ 2.4")
                            # simplest: do nothing special
                            # if you want to probe 0 explicitly (in same units as f*):
                            # if 0.0 not in scanned_points and 0.0 not in pending_sequence1:
                            #     pending_sequence1.insert(0, 0.0)

                        # If we decided on a refined center, build new scan and replace sequence
                        if new_center is not None:
                            print(f"Changing center to {new_center}")
                            new_sequence1 = self.make_scan_list(
                                center=new_center,
                                half_range_kHz=self.freq_scan_half_range_kHz,
                                min_step_kHz=self.freq_scan_step_size_kHz,
                                mode="pair",
                            )
                            pending_sequence1 = list(new_sequence1)
                            did_refine = True
                        else:
                            if all(r > 0.8 for r in (R0, R1, R2, R3)):
                                print("Resonance not in this range; seaching broader range")
                                new_sequence1 = self.make_scan_list(
                                    center=new_center,
                                    half_range_kHz=2*self.freq_scan_half_range_kHz,
                                    min_step_kHz=self.freq_scan_step_size_kHz,
                                    mode="pair",
                                )




            #### for fitting
            if self.enable_fitting and not self.scan_type.startswith("Ramsey"):
                BothSPCMs_RO1 = self.get_dataset("BothSPCMs_RO1")
                BothSPCMs_RO2 = self.get_dataset("BothSPCMs_RO2")
                retention_array = self.get_retention(BothSPCMs_RO1, BothSPCMs_RO2, self.n_measurements, len(self.scan_sequence1),
                                   self.single_atom_threshold * self.t_SPCM_first_shot)

                # print("inside fitting: ", retention_array)
                t = self.scan_sequence1
                y = retention_array

                p, p_err = self.which_fit_model.fit(t, y, evaluate_function=False, initialise= self.initialise)
                print(p)

                ### saving fitted parameters
                for var, val in p.items():
                    ds_name = f"fit_parameter_{var}"
                    self.set_dataset(ds_name, round(float(val),7), broadcast=True, persist=True)

                ### saving fitted parameters
                for var, val in p_err.items():
                    ds_name = f"fit_parameter_{var}_err"
                    self.set_dataset(ds_name, round(float(val),7), broadcast=True, persist=True)

                ### depending on scan_type, the fitted parameter is different
                if self.scan_type.startswith("Freq"):
                    if self.health_check_general(fit_check = True):
                        print("optimization success - dataset updating to fitted value")

                        self.set_dataset(scan_dict[self.scan_type]["center"], round(float(p["f0"]), -3), broadcast=True, persist=True)
                        print(scan_dict[self.scan_type]["center"]," updated to ", getattr(self, scan_dict[self.scan_type]["center"]))

                    else:
                        print("optimization failed - dataset not updated")

                elif self.scan_type.startswith("Time"):
                    ### updating dataset with fitted parameters
                    if self.health_check_general(fit_check=True):
                        print("optimization success - dataset updating to fitted value")
                        self.set_dataset(scan_dict[self.scan_type]["pi_pulse"], round(float(p["t_pi"]), 7), broadcast=True,
                                         persist=True)
                        print(scan_dict[self.scan_type]["pi_pulse"], " updated to ",
                              getattr(self, scan_dict[self.scan_type]["pi_pulse"]))

                    else:
                        print("optimization failed - dataset not updated")

                # write and overwrite the file here so we can quit the experiment early without losing data
                self.write_results({'name': self.experiment_name[:-11] + "_scan_over_" + self.scan_var_filesuffix})

                iteration += 1


        print("****************    General Variable Scan DONE   *****************")

    ################################### Health Check Functions ########################################

    def health_check_general(self, fit_check = False):
        """
        Perform a single-point health check on the microwave or time-domain scan.

        This method re-runs the experiment at the current resonance (or π-pulse)
        condition and evaluates whether the resulting fidelity meets the required
        `target_fidelity`. Depending on the selected scan type (frequency or time)
        and fit model (dip/peak), it computes the retention-based fidelity and
        updates the appropriate health-check dataset.

        Returns
        -------
        True  : if the health check passes (fidelity ≥ target_fidelity)
        False : if it fails and optimization is required
        """

        self.initialize_hardware()
        self.reset_datasets()

        if fit_check:
            if self.scan_type.startswith("Freq"):
                print("Original f0 value: ", getattr(self, scan_dict[self.scan_type]["center"]))
                fit_f0 = self.get_dataset("fit_parameter_f0")
                setattr(self, scan_dict[self.scan_type]["center"], round(float(fit_f0), -3))
                print("Fit check with fitted f0: ",  getattr(self, scan_dict[self.scan_type]["center"]))
            elif self.scan_type.startswith("Time"):
                print("Original t_microwave_00_pulse value: ", getattr(self, scan_dict[self.scan_type]["pi_pulse"]))
                fit_pi_pulse = self.get_dataset("fit_parameter_t_pi")
                setattr(self, scan_dict[self.scan_type]["pi_pulse"], round(float(fit_pi_pulse), 7))
                print("Fit check with fitted pi pulse: ", getattr(self, scan_dict[self.scan_type]["pi_pulse"]))
        else:
            # override specific variables. this will apply to the entire scan, so it is outside the loops
            for variable, value in self.override_ExperimentVariables_dict.items():
                setattr(self, variable, value)

            if self.scan_type.startswith("Freq"):
                print("Health Check with original f0: ",  getattr(self, scan_dict[self.scan_type]["center"]))
            elif self.scan_type.startswith("Time"):
                print("Health Check with original pi pulse: ", getattr(self, scan_dict[self.scan_type]["pi_pulse"]))


        # todo: make sure this does not interrupt the current dataset when implementing in GVS
        # override items:
        for var, val in scan_dict[self.scan_type]["override_items"].items():
            self.override_ExperimentVariables_dict[var] = getattr(self, val)

        if self.scan_type.startswith("Freq"):
            setattr(self, scan_dict[self.scan_type]["scan_variable1_name"], getattr(self, scan_dict[self.scan_type]["center"]))
        elif self.scan_type.startswith("Time"):
            setattr(self, scan_dict[self.scan_type]["scan_variable1_name"], getattr(self, scan_dict[self.scan_type]["pi_pulse"]))

        ### setting experiment function
        if self.which_node == 'bob':
            self.experiment_name = scan_dict[self.scan_type]["experiment_name_bob"]
        elif self.which_node == 'alice':
            self.experiment_name = scan_dict[self.scan_type]["experiment_name_alice"]

        ### setting target retention to calculate fidelity
        if scan_dict[self.scan_type]["fit_model"] in ["resonance_dip", "rabi_flop"]:
            target_retention_0 = True
        elif scan_dict[self.scan_type]["fit_model"] in ["resonance_peak", "rabi_flop_reversed"]:
            target_retention_0 = False

        self.initialise = scan_dict[self.scan_type]["initialise"]

        self.experiment_function = lambda: eval(self.experiment_name)(self)
        self.experiment_function()


        #todo: if this happens within some experiment, make sure it does not disturb the current dataset.
        BothSPCMs_RO1 = self.get_dataset("BothSPCMs_RO1")
        BothSPCMs_RO2 = self.get_dataset("BothSPCMs_RO2")

        # print("BothSPCMs_RO1", BothSPCMs_RO1)
        retention_array, loading_rate_array, n_atoms_loaded_array = self.get_loading_and_retention(BothSPCMs_RO1,
                                                                                                   BothSPCMs_RO2,
                                                                                                   self.n_measurements,
                                                                                                   int((len(BothSPCMs_RO1)-1)/(self.n_measurements)),
                                                                                                   self.single_atom_threshold * self.t_SPCM_first_shot)
        # print("retention_array", retention_array)
        # todo: break; if loading_rate too low

        if target_retention_0:
            fidelity = 1.0 - retention_array[-1]
        else:
            fidelity = retention_array[-1]

        #####update
        ##### update health_check_dataset only with FREQ scans;
        if self.scan_type.startswith("Freq"):
            self.set_dataset(scan_dict[self.scan_type]["health_check_dataset_name"], round(float(fidelity), 2), broadcast=True, persist=True)

        if fidelity >= self.target_fidelity:
            print("Health Check - passed with fidelity: ", fidelity)

            #### update health_check dataset with TIME scanse only when it passed the test.
            if self.scan_type.startswith("Time"):
                self.set_dataset(scan_dict[self.scan_type]["health_check_dataset_name"], round(float(fidelity), 2),
                                 broadcast=True, persist=True)
            return True
        else:
            return False

    def submit_opt_exp_general(self, override_arguments = None):  # todo: account for changes
        """
        Schedules an optimization experiment when a health check fails.

        Constructs a new `expid` dictionary for the optimization scan and injects
        the current scan configuration, override variables, and user-defined scan
        options. It then submits this optimization run to the ARTIQ scheduler with
        high priority (99), ensuring it executes before any follow-up scans.

        Parameters
        ----------
        override_arguments : dict, optional
            Key-value pairs to overwrite default optimization arguments.
            Example: {"enable_fitting": False}
        """
        print("submitting another experiment")
        ## 99 seems to be the highest priority that can be set to.

        # todo: make a default expid and overwrite it just a few things.
        default_expid = {
            "log_level": 30, #todo: check which level this is - debug? info? or else?
            "file": "qn_artiq_routines\\MicrowaveScanOptimizer.py",
            "class_name": "MicrowaveScanOptimizer",
            "arguments": {
                "run_health_check_and_optimize": False,
                "target_fidelity": 0.80,

                "n_measurements": self.n_measurements,
                "override_ExperimentVariables": "{'dummy_variable': 4}",
                "enable_fitting": True,
                "enable_geometric_frequency_scan": True,

                # Scan type toggles
                "Frequency_00_Scan": False,
                "Frequency_01_Scan": False,
                "Frequency_11_Scan": False,
                "Frequency_m10_Scan": False,
                "Time_00_Scan": False,
                "Time_01_Scan": False,
                "Time_11_Scan": False,
                "Time_m10_Scan": False,
                "Ramsey_00_Scan": False,
                "Ramsey_01_Scan": False,
                "Ramsey_11_Scan": False,

                # [Full Scan] Frequency scan parameters (kHz)
                "freq_scan_range_left_kHz": 100.0,
                "freq_scan_range_right_kHz": 101.0,
                "freq_scan_step_size_kHz": 20.0,

                # [Faster Scan] Frequency scan parameters - centered at resonance (kHz)
                "shrink_factor": 2.5,
                "freq_scan_half_range_kHz": 100.0,
                "freq_scan_min_step_size_kHz": 10.0,

                # # Time scan parameters
                "time_scan_sequence": 'np.arange(0,10,1)*us'

            },
            "repo_rev": "N/A",
        }

        new_expid = copy.deepcopy(default_expid)

        ###todo: update_default_expid_from_self - think about what should be updated based on initial;;

        for key in scan_options:
            if getattr(self, key, False):     # if self.key does not exist, it returns False
                new_expid["arguments"][key] = True
                break  # stop once the first TRUE is found

        ### keeping set override_ExperimentVariables
        new_expid["arguments"]["override_ExperimentVariables"] = self.override_ExperimentVariables


        if override_arguments is not None:
            for key, value in override_arguments.items():
                if key in new_expid["arguments"]:
                    print(f"Overriding {key}: {new_expid['arguments'][key]} to {value}")
                    new_expid["arguments"][key] = value

                else:
                    raise KeyError(f"Invalid override key: '{key}'. ")

        self.scheduler.submit(pipeline_name="main", expid=new_expid, priority=99, due_date=None, flush=False)

    def make_scan_list(self, center, half_range_kHz, min_step_kHz=10.0, method="center_geometric", mode="pair"):
        """
        Create a geometric scan list around a center frequency.

        * Supports two geometric scan strategies:
            - "center_geometric": symmetric shrinking steps toward the center
            - "quarter_geometric": concentrated sampling around ±half-range/2 clusters

        * Two sequencing modes define point ordering:
            - "sequential": sorted dense sampling from edges toward the center
            - "pair": symmetric (-,+) probing at each geometric shell

        ** Uses a shrink ratio r = (1 - 1/self.shrink_factor), with self.shrink_factor > 1.

        Modes
        -----
        mode = "sequential":
            offsets: [x, x*r, x*r^2, ... >= min_step_kHz]
            final list (sorted):
                [-x, -x*r, -x*r^2, ..., 0, ..., x*r^2, x*r, x]

        mode = "pair":
            offsets in order:
                [-x, +x, -x*r, +x*r, -x*r^2, +x*r^2, ..., 0]

        """
        assert self.shrink_factor > 1, "shrink_factor must be > 1"
        assert mode in ("sequential", "pair"), f"Unknown mode: {mode}"

        x = float(half_range_kHz)
        r = (1 - 1 / self.shrink_factor)

        # if method == "center_geometric":
        #     # Build offsets depending on mode
        #     if mode == "sequential":
        #         offsets = []
        #         val = x
        #         while val >= min_step_kHz:
        #             offsets.append(val)
        #             val *= r
        #
        #         symmetric_offsets = sorted([-v for v in offsets] + [0.0] + offsets)
        #         values = [center + v * kHz for v in symmetric_offsets]
        #
        #     elif mode == "pair":
        #         offsets = []
        #         val = x
        #         while val >= min_step_kHz:
        #             offsets.append(-val)
        #             offsets.append(+val)
        #             val *= r
        #
        #         offsets.append(0.0)
        #         values = [center + v * kHz for v in offsets]


        if method == "center_geometric":
            if mode == "sequential":
                offsets = []
                val = x
                prev = None
                while val >= min_step_kHz and (prev is None or abs(prev - val) >= min_step_kHz):
                    offsets.append(val)
                    prev = val
                    val *= r

                symmetric_offsets = sorted([-v for v in offsets] + [0.0] + offsets)
                values = [center + v * kHz for v in symmetric_offsets]

            elif mode == "pair":
                offsets = []
                val = x
                prev = None
                while val >= min_step_kHz and (prev is None or abs(prev - val) >= min_step_kHz):
                    offsets.append(-val)
                    offsets.append(+val)
                    prev = val
                    val *= r

                offsets.append(0.0)
                values = [center + v * kHz for v in offsets]


        elif method == "quarter_geometric":
            ### quarter span (in kHz) – clusters will be centered at ±quarter
            quarter = x / 2.0

            # Local "center_geometric" grid around 0 in [-quarter, +quarter]
            local_pos = []
            val_local = quarter
            while val_local >= min_step_kHz:
                local_pos.append(val_local)
                val_local *= r

            local_offsets = sorted([-v for v in local_pos] + [0.0] + local_pos)

            if mode == "sequential":
                raw_offsets = []

                # Left cluster: centered at -quarter
                for o in local_offsets:
                    raw_offsets.append(-quarter + o)

                # Right cluster: centered at +quarter
                for o in local_offsets:
                    raw_offsets.append(+quarter + o)

                # Always include global center
                raw_offsets.append(0.0)

                # Clip and deduplicate
                offsets = sorted({o for o in raw_offsets if -x <= o <= x})
                values = [center + v * kHz for v in offsets]

            elif mode == "pair":
                raw_offsets = []
                raw_offsets.append(-quarter)
                raw_offsets.append(+quarter)

                for o in local_offsets:
                    if o == 0.0:
                        continue
                    candidates = [
                        -quarter + o,
                        -quarter - o,
                        +quarter + o,
                        +quarter - o,
                    ]
                    for c in candidates:
                        if -x <= c <= x:
                            raw_offsets.append(c)

                if 0.0 not in raw_offsets:
                    raw_offsets.append(0.0)

                offsets = raw_offsets
                values = [center + v * kHz for v in offsets]


        return np.array(values)

    def get_scan_type(self):
        """
        Determine which scan type (e.g., Frequency_00_Scan) is currently enabled.

        Ensures that exactly one scan flag is True among `scan_options`. Prevents
        ambiguous or conflicting scan configuration before running the experiment.
        """

        # Collect all scan flags that are True
        enabled = [name for name in scan_options if getattr(self, name, False)]

        if len(enabled) == 0:
            raise ValueError("No scan type selected. At least one scan flag must be True.")

        if len(enabled) > 1:
            raise ValueError(
                f"Multiple scan types selected ({enabled}). Only one can be True."
            )

        # Exactly one scan type selected
        return enabled[0]

    def update_default_expid_from_self(self, expid):
        """
        Update expid['arguments'] from attributes on self.
        Assumes all keys already exist and are set on self.
        """
        args = expid['arguments']
        for key in list(args.keys()):
            if hasattr(self, key):
                args[key] = getattr(self, key)

        return expid

    def get_loading_and_retention(self, photocounts, photocounts2, measurements, iterations, cutoff):
        """
        Compute loading rate, retention, and number of atoms loaded per iteration.

        Returns:
        * retention_array[iterations]
        * loading_rate_array[iterations]
        * n_atoms_loaded_array[iterations]
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
        Compute only the retention fraction per iteration.

        - Same logic as get_loading_and_retention(), but returns only
          retention_array[iterations].
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

