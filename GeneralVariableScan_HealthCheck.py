"""
For scanning a defined experiment function over an ExperimentVariable

In many cases, we want the flexibility to be able to scan an experiment over a wide range of parameters--
in principle, over any of the defined ExperimentVariables. Including all of these as ARTIQ ScanVariables
in the GUI would be cumbersome for one experiment, not to mention including this for several experiments.
This code allows the user to scan up in up to 2 dimensions by supplying ExperimentVariables by name,
corresponding python sequences defining the scan steps, and an experiment function defined in
utilities/experiment_functions.py. As some of these variables pertain to hardware settings, such as DDS power,
it is necessary in general to re-initialize hardware at each scan step. We accomplish this by calling
base.build before each call of the experiment function, i.e., at the start of each scan step.
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

from GeneralVariableScan_Microwaves import scan_dict


class GeneralVariableScan_HealthCheck(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        ## health check
        self.setattr_argument('run_health_check_and_schedule', BooleanValue(default=False), "Health Check")
        self.setattr_argument("target_fidelity", NumberValue(0.80, ndecimals=2, step=1), "Health Check")
        self.setattr_argument("health_check_with_n_measurements", NumberValue(100, ndecimals=0, step=1), "Health Check")

        self.setattr_argument('health_check_every_n_ite', BooleanValue(default=False), "Health Check")
        self.setattr_argument('health_check_every_delta_t_hours', BooleanValue(default=False), "Health Check")

        self.setattr_argument("every_n_ite", NumberValue(10, ndecimals=0, step=1), "Health Check")
        self.setattr_argument("every_delta_t_hours", NumberValue(0.80, ndecimals=1, step=1), "Health Check")

        #todo: health check microwaves freq scans - multiple - ["","",""] as a list? or boolean?


        # the number of measurements to be made for a certain setting of the
        # experiment parameters
        self.setattr_argument("n_measurements", NumberValue(100, ndecimals=0, step=1), "General Variable Scan")
        self.setattr_argument('scan_variable1_name', StringValue('t_blowaway'), "General Variable Scan")
        self.setattr_argument("scan_sequence1", StringValue(
            'np.array([0.000,0.005,0.02,0.05])*ms'), "General Variable Scan")

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

        self.override_ExperimentVariables_dict = eval(self.override_ExperimentVariables)
        assert type(self.override_ExperimentVariables_dict) == dict, \
            "override_ExperimentVariables should be a python dictionary"

        for variable, value in self.override_ExperimentVariables_dict.items():
            assert hasattr(self, variable), (f"There is no ExperimentVariable " + variable +
                                                    ". Did you mistype it?")

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

                self.initialize_dependent_variables()
                self.initialize_hardware()
                self.reset_datasets()

                # the measurement loop.
                self.experiment_function()


                # write and overwrite the file here so we can quit the experiment early without losing data
                self.write_results({'name': self.experiment_name[:-11] + "_scan_over_" + self.scan_var_filesuffix})

                iteration += 1

                ## Run Health Check
                if self.run_health_check_and_schedule:
                    if self.health_check_every_n_ite and (iteration % self.every_n_ite) == 0:
                        print(f"running health check after iteration # {iteration-1}")

                        ### run health check
                        failed_scans = self.health_check_microwave_freqs()

                        if failed_scans:
                            print("These scans need re-optimisation:", failed_scans)
                            ### if health check fail, i) schedule optimization

                            ### if health check fail, ii) schedule resuming experiment

                        else:
                            print("All microwave health checks passed.")

                            ### update everything back to previous setting


                            ### Override variables again to avoid conflict
                            # override specific variables. this will apply to the entire scan, so it is outside the loops
                            for variable, value in self.override_ExperimentVariables_dict.items():
                                setattr(self, variable, value)

                            ### resume the current scan



        print("****************    General Variable Scan DONE   *****************")


    def health_check_microwave_freqs(self):

        prev_exp_fun = self.experiment_name


        scan_options = ["Frequency_00_Scan"]
        # scan_options = ["Frequency_00_Scan", "Frequency_01_Scan", "Frequency_11_Scan"]

        for i, scan_type in enumerate(scan_options):

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
            print("Running experiment...")


            self.experiment_function = lambda: eval(self.experiment_name)(self)
            self.experiment_function()
            print("Experiment finished.")


            # todo: make sure it does not disturb the current dataset.... - does not affect the actual dataset
            BothSPCMs_RO1 = self.get_dataset("BothSPCMs_RO1")
            BothSPCMs_RO2 = self.get_dataset("BothSPCMs_RO2")

            # print("BothSPCMs_RO1", BothSPCMs_RO1)
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
                # self.n_measurements = record_prev_val
                # self.set_dataset("n_measurements", record_prev_val, persist=True)

            else:
                print(f"Health Check {scan_type} - failed with fidelity: {fidelity}")
                return scan_options[i:]

        ### if passed the test: update everything back to previous setting
        self.experiment_function = lambda: eval(prev_exp_fun)(self)

        return []



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