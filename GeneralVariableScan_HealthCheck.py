"""
GeneralVariableScan_HealthCheck

This class runs a GVS with chosen experiment function, periodically performs
single-point microwave health checks, and if a check fails, it automatically schedules a dedicated
optimization scan followed by a resumed health-check scan that continues the original sequence from
where it left off, while keeping all necessary ExperimentVariables and overrides consistent across experiments.

Note - cannot change n_measurements just for health_check because it is a non-persistent variable.
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

from GeneralVariableScan_Microwaves import scan_dict, scan_options

import copy

class GeneralVariableScan_HealthCheck(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        ## For Bookeeping purpose
        self.setattr_argument("parent_rid", NumberValue(0, ndecimals=0, step=1), "** For bookeeping purpose - always set this to 0 **")

        ## health check
        self.setattr_argument('run_health_check_and_schedule', BooleanValue(default=True), "Health Check")
        self.setattr_argument("target_fidelity", NumberValue(0.80, ndecimals=2, step=1), "Health Check")
        self.setattr_argument('override_arguments_for_scheduling_optimization', StringValue("{'n_measurements':100}"),
                              "Health Check")
        # self.setattr_argument("health_check_with_n_measurements", NumberValue(100, ndecimals=0, step=1), "Health Check")

        self.setattr_argument('health_check_every_n_ite', BooleanValue(default=True), "Health Check")
        self.setattr_argument('health_check_every_delta_t_hours', BooleanValue(default=False), "Health Check - not yet implemented")

        self.setattr_argument("every_n_ite", NumberValue(5, ndecimals=0, step=1), "Health Check")
        self.setattr_argument("every_delta_t_hours", NumberValue(0.80, ndecimals=1, step=1), "Health Check - not yet implemented")

        self.setattr_argument("Frequency_00_Scan", BooleanValue(default=False),"Health Check - Select to include in health check")
        self.setattr_argument("Frequency_01_Scan", BooleanValue(default=False),"Health Check - Select to include in health check")
        self.setattr_argument("Frequency_11_Scan", BooleanValue(default=False),"Health Check - Select to include in health check")
        self.setattr_argument("Frequency_m10_Scan", BooleanValue(default=False),"Health Check - Select to include in health check")


        # the number of measurements to be made for a certain setting of the
        # experiment parameters
        self.setattr_argument("n_measurements", NumberValue(100, ndecimals=0, step=1), "General Variable Scan")
        self.setattr_argument('scan_variable1_name', StringValue('dummy_variable'), "General Variable Scan")
        self.setattr_argument("scan_sequence1", StringValue(
            'np.arange(0,5,1)'), "General Variable Scan")

        # this variable is optional
        self.setattr_argument('scan_variable2_name', StringValue(''), "General Variable Scan")
        self.setattr_argument("scan_sequence2", StringValue(
            'np.linspace(-2,2,5)*V'), "General Variable Scan")

        # allows user to supply a dictionary of values to override. this is useful for when
        # you don't want to constantly go check ExperimentVariables to see if, e.g. blowaway_light_off
        # is False. You can just set it here to guarantee the behavior you want, without changing
        # the value stored in the dataset, so subsequent experiments will be unaffected. This leads
        # to fewer errors overall.
        self.setattr_argument('override_ExperimentVariables', StringValue("{'dummy_variable':4}"), "General Variable Scan")

        experiment_function_names_list = [x for x in dir(exp_functions)
            if ('__' not in x and str(type(getattr(exp_functions,x)))=="<class 'function'>"
                and 'experiment' in x)]

        # a function that take no arguments that gets imported and run
        self.setattr_argument('experiment_function',
                              EnumerationValue(experiment_function_names_list), "General Variable Scan")

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

        ### If this is the very first run (user-submitted), parent_rid will be 0
        if int(self.parent_rid) == 0:
            self.parent_rid = int(self.scheduler.rid)
            print("[Parent initialized] This run is the parent. Setting parent_rid =", self.parent_rid)
        else:
            print("[Child scan] Inheriting parent_rid =", self.parent_rid)

        self.set_dataset("parent_rid", self.parent_rid, broadcast=True, persist=True)

        self.scan_variable1 = str(self.scan_variable1_name)
        self.scan_variable2 = str(self.scan_variable2_name)

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

        # handle some special scan cases
        if self.scan_variable1 == 'f_Rigol_modulation':

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
        if self.scan_variable2 == 'f_Rigol_modulation':
            self.print_async("why would you use scan_variable2 for this?")
            raise

        scan_vars = [self.scan_variable1_name, self.scan_variable2_name]
        scan_vars = [x for x in scan_vars if x != '']
        self.scan_var_labels = ','.join(scan_vars)
        self.scan_var_filesuffix = '_and_'.join(scan_vars)

        ### override_ExperimentVariables_dict for GVS
        self.override_ExperimentVariables_dict = eval(self.override_ExperimentVariables)
        assert type(self.override_ExperimentVariables_dict) == dict, \
            "override_ExperimentVariables should be a python dictionary"

        for variable, value in self.override_ExperimentVariables_dict.items():
            assert hasattr(self, variable), (f"There is no ExperimentVariable " + variable +
                                                    ". Did you mistype it?")


        ### override_arguments_for_scheduling_optimization for scheduling optimization experiment
        self.override_arguments_for_scheduling_optimization_dict = eval(self.override_arguments_for_scheduling_optimization)
        assert type(self.override_arguments_for_scheduling_optimization_dict) == dict, \
            "override_arguments_for_scheduling_optimization_dict should be a python dictionary"

        # for variable, value in self.override_arguments_for_scheduling_optimization_dict.items():
        #     assert hasattr(self, variable), (f"There is no ExperimentVariable " + variable +
        #                                             ". Did you mistype it?")


        try:
            self.experiment_name = self.experiment_function
            self.experiment_function = lambda :eval(self.experiment_name)(self)
        except NameError as e:
            print(f"The function {self.experiment_name} is not defined. Did you forget to import it?")
            raise

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

        #print(status_dict)

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
        Main scan loop with embedded health checks.

        - Optionally rebuilds/repairs base experiment if datasets may be stale.
        - Initializes datasets and applies override_ExperimentVariables.
        - Loops over scan_sequence1 and scan_sequence2:
            * Sets scan variables.
            * Re-prepares dependent variables and hardware.
            * Runs the selected experiment_function.
            * Writes scan results to file.
        - Every `every_n_ite` iterations (if enabled), runs a microwave
          health check:
            * If all checks pass → resumes scan.
            * If any checks fail → schedules optimization scan(s) and a
              "resume" GeneralVariableScan_HealthCheck, then terminates
              the current run and pauses the scheduler.
        """

        if self.needs_fresh_build:
            self.base.build()
            self.base.prepare()

        self.initialize_datasets()

        # scan in up to 2 dimensions. for each setting of the parameters, run experiment_function n_measurement times
        iteration = 0
        self.set_dataset("iteration", iteration, broadcast=True)
        # self.set_dataset("n_measurements", self.n_measurements, persist=True)

        # override specific variables. this will apply to the entire scan, so it is outside the loops
        for variable, value in self.override_ExperimentVariables_dict.items():
            setattr(self, variable, value)

        # self.warm_up()

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

                self.initialize_dependent_variables() ## base.prepare()
                self.initialize_hardware()
                self.reset_datasets()

                # print("experiment_function running with n_measurements :", self.n_measurements)

                # the measurement loop.
                self.experiment_function()


                # write and overwrite the file here so we can quit the experiment early without losing data
                self.write_results({'name': self.experiment_name[:-11] + "_scan_over_" + self.scan_var_filesuffix})

                iteration += 1

                ## Run Health Check
                if self.run_health_check_and_schedule:
                    if self.health_check_every_n_ite and (iteration % self.every_n_ite) == 0:
                        print(f"running health check after iteration # {iteration-1}")


                        ####todo:  error if nothing checked;
                        ### run health check
                        failed_scans = self.health_check_microwave_freqs()

                        # write and overwrite the health check results
                        self.write_results(
                            {'name': self.experiment_name[:-11] + "_scan_over_" + self.scan_var_filesuffix})

                        if failed_scans:
                            print("These scans need re-optimisation:", failed_scans)
                            ###todo: if health check fail, i) schedule optimization
                            for scan_name in failed_scans:
                                print("Scheduling experiment to optimize: ", scan_name)
                                override_args = copy.deepcopy(self.override_arguments_for_scheduling_optimization_dict)
                                override_args[scan_name] = True
                                self.submit_optimization_scans(override_arguments = override_args)
                                self.override_arguments_for_scheduling_optimization_dict[scan_name] = False

                            ###todo: if health check fail, ii) schedule resuming experiment
                            self.submit_resume_scan_after_optimization(current_iteration = iteration-1)

                            ### after scheduling scans above, it terminates the current experiment
                            self.scheduler.request_termination(self.scheduler.rid)

                        else:
                            print("All microwave health checks passed. Resuming next scan")

                            ### update everything back to previous setting - done inside healthcheck if passed

                            ### Override variables again to avoid conflict
                            # override specific variables. this will apply to the entire scan, so it is outside the loops
                            for variable, value in self.override_ExperimentVariables_dict.items():
                                setattr(self, variable, value)

                            ### resume the current scan
                        ### terminating the experiment if health_check failed.
                        self.core.comm.close()  # placing the hardware in a safe state and disconnecting it from the core device
                        self.scheduler.pause()



        print("****************    General Variable Scan DONE   *****************")

    ################################### Health Check Functions ########################################
    #### slightly different from health check functions in GVS_Microwaves.

    def health_check_microwave_freqs(self):
        """
        Run single-point microwave health checks for selected scans.

        - Looks at which Frequency_*_Scan flags are True (in scan_options).
        - For each enabled scan_type:
            * Sets up override_ExperimentVariables for that scan.
            * Chooses the correct experiment_function (alice/bob).
            * Runs one health-check experiment at the nominal center frequency.
            * Computes retention using get_loading_and_retention().
            * Converts retention to a fidelity metric and writes it to a
              dedicated health-check dataset.
            * If fidelity < target_fidelity:
                - returns the remaining enabled scan types (this and later ones)
                  as `failed_scans`.

        - If all scans pass, returns an empty list.

        """
        prev_exp_fun = self.experiment_name

        enabled_scan_options = [name for name in scan_options if getattr(self, name, False)]
        print("enabled_scan_options ", enabled_scan_options)

        for i, scan_type in enumerate(enabled_scan_options):

            print("Health Check - ", scan_type)
            self.scan_type = scan_type

            self.initialize_hardware()
            self.reset_datasets()

            # override specific variables. this will apply to the entire scan, so it is outside the loops
            for variable, value in self.override_ExperimentVariables_dict.items():
                setattr(self, variable, value)

            print("Health Check with original f0: ", getattr(self, scan_dict[self.scan_type]["center"]))

            # todo: make sure this does not interrupt the current dataset when implementing in GVS
            # override items:
            for var, val in scan_dict[self.scan_type]["override_items"].items():
                self.override_ExperimentVariables_dict[var] = getattr(self, val)


            setattr(self, scan_dict[self.scan_type]["scan_variable1_name"],
                        getattr(self, scan_dict[self.scan_type]["center"]))


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


            # record_prev_val = self.n_measurements
            # print(self.n_measurements)
            # self.n_measurements = self.health_check_with_n_measurements
            # self.set_dataset("n_measurements", self.health_check_with_n_measurements, persist=True)

            # print(self.n_measurements)
            # print("self.target_fidelity in health_check within GVShealthcheck", self.target_fidelity)
            print("Running experiment...")
            self.in_health_check = True
            # print("self.in_health_check - before exp ",self.in_health_check)

            self.experiment_function = lambda: eval(self.experiment_name)(self)
            self.experiment_function()

            ### update the experiment function back to previous setting
            self.experiment_name = prev_exp_fun
            self.experiment_function = lambda: eval(prev_exp_fun)(self)


            # todo: make sure it does not disturb the current dataset.... - does not affect the actual dataset
            BothSPCMs_RO1 = self.get_dataset("BothSPCMs_RO1_in_health_check")
            BothSPCMs_RO2 = self.get_dataset("BothSPCMs_RO2_in_health_check")

            retention_array, loading_rate_array, n_atoms_loaded_array = self.get_loading_and_retention(BothSPCMs_RO1,
                                                                                                       BothSPCMs_RO2,
                                                                                                       self.n_measurements,
                                                                                                       int((len(BothSPCMs_RO1) - 1) / (self.n_measurements)),
                                                                                                       self.single_atom_threshold * self.t_SPCM_first_shot)
            # print("retention_array", retention_array)
            # todo: break; if loading_rate too low

            if target_retention_0:
                fidelity = 1.0 - retention_array[-1]
            else:
                fidelity = retention_array[-1]

            #####update
            self.set_dataset(scan_dict[self.scan_type]["health_check_dataset_name"], round(float(fidelity), 2),
                                 broadcast=True, persist=True)

            if fidelity >= self.target_fidelity:
                print(f"Health Check {scan_type} - passed with fidelity: {fidelity}")

            else:
                print(f"Health Check {scan_type} - failed with fidelity: {fidelity}")

                return enabled_scan_options[i:]

        self.in_health_check = False

        return []

    def submit_optimization_scans(self, override_arguments = None):
        """
        Schedule a GeneralVariableScan_Microwaves optimization experiment.

        - Builds a default expid for GeneralVariableScan_Microwaves with:
            * run_health_check_and_optimize = False
            * enable_faster_frequency_scan = True
            * all Frequency_*_Scan flags initially False
            * default frequency scan parameters.
        - Replaces override_ExperimentVariables with the current
          override_ExperimentVariables_dict (stringified).
        - Applies `override_arguments` (e.g. enabling a specific Frequency_*_Scan,
          changing n_measurements, etc.), checking keys against default arguments.
        - Submits the optimization job to the scheduler with priority 99.
        """
        print("submitting optimization experiment...")
        ## 99 is the highest priority that can be set to.
        ## making a default expid and overwriting it just a few things.

        # todo: make a default expid and overwrite it just a few things.
        default_expid = {
            "log_level": 30,
            "file": "qn_artiq_routines\\GeneralVariableScan_Microwaves.py",
            "class_name": "GeneralVariableScan_Microwaves",
            "arguments": {
                "run_health_check_and_optimize": False,
                "target_fidelity": 0.8,
                "n_measurements": 100,
                "override_ExperimentVariables": "{'dummy_variable':4}",

                # scan control
                "enable_faster_frequency_scan": True,
                "enable_fitting": True,

                # which scans to run
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

                # frequency scan parameters
                "freq_scan_range_left_kHz": 100.0,
                "freq_scan_range_right_kHz": 100.0,
                "freq_scan_step_size_kHz": 10.0,
                "shrink_factor": 2.5,
                "freq_scan_half_range_kHz": 100.0,
                "freq_scan_min_step_size_kHz": 10.0,

                # time scan parameters
                "time_scan_sequence": "np.arange(0,10,.5)*us",
            },
            "repo_rev": "N/A",
        }

        new_expid = copy.deepcopy(default_expid)

        ### keeping set override_ExperimentVariables
        # new_expid["arguments"]["override_ExperimentVariables"] = self.override_ExperimentVariables

        override_ExperimentVariables_dict_str = repr(self.override_ExperimentVariables_dict)
        # print("override_ExperimentVariables_dict_str ", override_ExperimentVariables_dict_str)

        new_expid["arguments"]["override_ExperimentVariables"] = override_ExperimentVariables_dict_str
        new_expid["arguments"]["target_fidelity"] = self.target_fidelity

        print("override_arguments in sumbmit_optimization_scans: ", override_arguments)

        if override_arguments is not None:
            for key, value in override_arguments.items():
                if key in new_expid["arguments"]:
                    print(f"Overriding {key}: {new_expid['arguments'][key]} to {value}")
                    new_expid["arguments"][key] = value

                else:
                    raise KeyError(f"Invalid override key: '{key}'. ")



        self.scheduler.submit(pipeline_name="main", expid=new_expid, priority=99, due_date=None, flush=False)

    def submit_resume_scan_after_optimization(self, current_iteration, override_arguments = None):  # todo: account for changes
        """
        Schedule a follow-up GeneralVariableScan_HealthCheck to resume a scan.

        - Builds a new expid for GeneralVariableScan_HealthCheck with default
          health-check and scan settings.
        - Copies relevant arguments from the current experiment via
          update_default_expid_from_self().
        - Constructs a new scan_sequence1 starting from `current_iteration`
          onward, so the resumed job continues where the previous one stopped.
        - Sets experiment_function to the current experiment_name (string).
        - Optionally applies additional `override_arguments`.
        - Serializes current override_ExperimentVariables_dict back into
          override_ExperimentVariables.
        - Submits the resumed scan to the scheduler with priority 98, so it
          runs after the optimization scans.
        """
        print("submitting another experiment that starts from iteration ", current_iteration)
        ## 99 seems to be the highest priority that can be set to.

        # todo: make a default expid and overwrite it just a few things.
        default_expid = {
            'log_level': 30,
            'file': 'qn_artiq_routines\\GeneralVariableScan_HealthCheck.py',
            'class_name': 'GeneralVariableScan_HealthCheck',
            'arguments': {
                'parent_rid': self.parent_rid,
                'run_health_check_and_schedule': True,
                'target_fidelity': 0.8,
                # overrides
                'override_arguments_for_scheduling_optimization': "{'enable_fitting': True}",
                'health_check_every_n_ite': True,
                'health_check_every_delta_t_hours': False,
                'every_n_ite': 2,
                'every_delta_t_hours': 0.8,

                # scan types
                'Frequency_00_Scan': False,
                'Frequency_01_Scan': False,
                'Frequency_11_Scan': False,
                'Frequency_m10_Scan': False,

                # measurement settings
                'n_measurements': 100,

                # scan variables
                'scan_variable1_name': 'dummy_variable',
                'scan_sequence1': 'np.arange(0,3.1,1)',
                'scan_variable2_name': '',
                'scan_sequence2': 'np.linspace(-2,2,5)*V',

                # overrides
                'override_ExperimentVariables': "{'dummy_variable': 4}",

                # experiment selection
                'experiment_function': 'atom_loading_2_experiment',
                'control_experiment': False,
            },
            'repo_rev': 'N/A'
        }

        new_expid = copy.deepcopy(default_expid)
        new_expid = self.update_default_expid_from_self(new_expid)

        ### New scan sequence - starting from the last iteration
        # todo: start from last iteration???? or something else????
        new_seq = self.scan_sequence1[current_iteration:]
        new_scan_sequence1 = f"np.array({new_seq.tolist()})"  # convert to numpy-based string expression
        # print("new_scan_sequence1: ", new_scan_sequence1)

        # write back to expid arguments
        new_expid["arguments"]["scan_sequence1"] = new_scan_sequence1

        # experiment function should be in string;
        new_expid["arguments"]["experiment_function"] = self.experiment_name

        if override_arguments is not None:
            for key, value in override_arguments.items():
                if key in new_expid["arguments"]:
                    print(f"Overriding {key}: {new_expid['arguments'][key]} to {value}")
                    new_expid["arguments"][key] = value

                else:
                    raise KeyError(f"Invalid override key: '{key}'. ")

        # print("sumbmit_resume - self.override_ExperimentVariables_dict: ", self.override_ExperimentVariables_dict)
        override_ExperimentVariables_dict_str = repr(self.override_ExperimentVariables_dict)
        new_expid["arguments"]["override_ExperimentVariables"] = override_ExperimentVariables_dict_str

        print("sumbmit_resume - override_ExperimentVariables_dict_str: ", override_ExperimentVariables_dict_str)

        child_rid = self.scheduler.submit(pipeline_name="main", expid=new_expid, priority=98, due_date=None, flush=False)

        self.set_dataset("child_rid", child_rid, broadcast=True, persist=True)
        print("Remaining scan resuming at RID: ", child_rid)

    def update_default_expid_from_self(self, expid):
        """
        Populate expid['arguments'] from attributes on this experiment.

        - For each key in expid['arguments']:
            * If this object (`self`) has an attribute with that name,
              copy its value into the argument dict.
        - Returns the modified expid.

        Used to carry over settings (n_measurements, scan flags, etc.) from
        the current HealthCheck experiment into newly scheduled jobs.
        """
        args = expid['arguments']
        for key in list(args.keys()):
            if hasattr(self, key):
                args[key] = getattr(self, key)

        return expid


    def get_loading_and_retention(self, photocounts, photocounts2, measurements, iterations, cutoff):
        """
        Compute loading rate, retention, and number of atoms loaded per iteration.

        - Splits photocounts into (iterations × measurements) blocks.
        - For each iteration:
            * Counts how many shots exceed `cutoff` in RO1 (atoms loaded).
            * Counts how many of those are still above cutoff in RO2
              (atoms retained).
            * Computes loading_rate = n_loaded / measurements.
            * Computes retention = n_retained / n_loaded (or 0 if no atoms).
        - Returns:
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