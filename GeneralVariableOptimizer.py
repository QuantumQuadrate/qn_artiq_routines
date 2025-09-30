"""
For optimizing an experiment function over any number of defined variables

Any experiment function defined in experiment_functions.py can be optimized
over any set of ExperimentVariables

Example strings for the variables_and_bounds argument:
    [('AZ_bottom_volts_PGC',-0.4*V,0.0*V,'diff'),('AZ_top_volts_PGC',-0.1*V,0.3*V,'diff'),
('AX_volts_PGC',-0.1*V,0.3*V,'diff'),('AY_volts_PGC',-0.1*V,0.3*V,'diff')]

    [('f_cooling_DP_RO',-2*MHz,2*MHz,'diff'),('p_cooling_DP_RO',0.75,1.0,'abs')]
    [('p_cooling_DP_PGC',0.1,1.0,'abs'),('f_cooling_DP_PGC',115*MHz,125*MHz,'abs'),
('t_PGC_after_loading',1*ms,40*ms,'abs')]

    [('set_point_PD1_AOM_A1',0.95,1.05,'perc'),('set_point_PD2_AOM_A2',0.95,1.05,'perc'),
('set_point_PD3_AOM_A3',0.95,1.05,'perc'),('set_point_PD4_AOM_A4',0.95,1.05,'perc'),
('set_point_PD5_AOM_A5',0.95,1.05,'perc'),('set_point_PD6_AOM_A6',0.95,1.05,'perc'),
('AZ_bottom_volts_MOT',-0.05*V,0.05*V,'diff'),('AZ_top_volts_MOT',-0.05*V,0.05*V,'diff'),
('AX_volts_MOT',-0.05*V,0.05*V,'diff'),('AY_volts_MOT',-0.05*V,0.05*V,'diff')]

See also GeneralVariableScan and AtomLoadingOptimizerMLOOP
"""
import logging

from artiq.experiment import *
import numpy as np
import scipy as sp
import logging

#Imports for M-LOOP
import mloop.interfaces as mli
import mloop.controllers as mlc
import mloop.visualizations as mlv

import sys, os
cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")
from utilities.BaseExperiment import BaseExperiment

# this is where your experiment function should live
from subroutines.experiment_functions import *
from subroutines.cost_functions import *
import subroutines.experiment_functions as exp_functions
import subroutines.cost_functions as cost_functions


class OptimizerVariable:
    """
    a container for the variables that we want to optimize.
    for streamlining setup
    """
    def __init__(self, var_and_bounds_tuple, experiment):
        self.name = var_and_bounds_tuple[0]
        assert hasattr(experiment, self.name), (f"There is no ExperimentVariable " + self.name +
                                                ". Did you mistype it?")
        self.is_differential = 1 if var_and_bounds_tuple[3] == 'diff' else 0
        self.is_percentage = 1 if var_and_bounds_tuple[3] == 'perc' else 0
        self.default_value = getattr(experiment, self.name)
        if self.is_percentage:
            self.min_bound = var_and_bounds_tuple[1] * self.default_value
            self.max_bound = var_and_bounds_tuple[2] * self.default_value
        else:
            self.min_bound = var_and_bounds_tuple[1] + self.is_differential*self.default_value
            self.max_bound = var_and_bounds_tuple[2] + self.is_differential*self.default_value


# Declare your custom class that inherits from the Interface class
class MLOOPInterface(mli.Interface):

    # Initialization of the interface, including this method is optional
    def __init__(self):
        # You must include the super command to call the parent class, Interface, constructor
        super(MLOOPInterface, self).__init__()

    # You must include the get_next_cost_dict method in your class
    # this method is called whenever M-LOOP wants to run an experiment
    def get_next_cost_dict(self, params_dict):  # we'll redefine this function later
        pass


class GeneralVariableOptimizer(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.setattr_argument("n_measurements", NumberValue(100, type='int', scale=1, ndecimals=0, step=1))
        self.setattr_argument("set_best_parameters_at_finish", BooleanValue(True))

        # [('experiment_var_name',minval,maxval,differential or absolute)]
        self.setattr_argument("variables_and_bounds", StringValue(
            "[('dummy_variable',0*ms,100*ms,'abs'),('dummy_variable',-0.5*V,0.5*V,'diff')]"))

        experiment_function_names_list = [x for x in dir(exp_functions)
            if ('__' not in x and str(type(getattr(exp_functions, x))) == "<class 'function'>"
                and 'experiment' in x)]

        cost_function_names_list = [x for x in dir(cost_functions)
            if ('__' not in x and str(type(getattr(cost_functions, x))) == "<class 'function'>"
                and 'cost' in x)]

        # a function that take no arguments that gets imported and run
        self.setattr_argument('experiment_function',
                              EnumerationValue(experiment_function_names_list, default='atom_loading_experiment'))

        # a function that takes a 1D sequence of photocounts data that gets imported and run
        self.setattr_argument('cost_function',
                              EnumerationValue(cost_function_names_list, default='atom_retention_cost'))

        # allows user to supply a dictionary of values to override. this is useful for when
        # you don't want to constantly go check ExperimentVariables to see if, e.g. blowaway_light_off
        # is False. You can just set it here to guarantee the behavior you want, without changing
        # the value stored in the dataset, so subsequent experiments will be unaffected. This leads
        # to fewer errors overall.
        self.setattr_argument('override_ExperimentVariables', StringValue("{'require_atom_loading_to_advance':False}"))

        group1 = "optimizer settings"
        self.setattr_argument("max_runs", NumberValue(70, type='int', scale=1, ndecimals=0, step=1), group1)
        self.setattr_argument("target_cost", NumberValue(-100, type='int', scale=1, ndecimals=0, step=1), group1)

        # todo: add keyword arguments for setting up M-LOOP
        #  there is a lot more available than what we've been using

        self.base.set_datasets_from_gui_args()
        logging.debug("build - done")

    def prepare(self):
        self.base.prepare()

        self.single_atom_SPCM0_RO1_threshold = self.single_atom_threshold*self.t_SPCM_first_shot
        self.single_atom_SPCM0_RO2_threshold = self.single_atom_threshold*self.t_SPCM_second_shot
        self.single_atom_RO1_threshold = self.single_atom_threshold*self.t_SPCM_first_shot
        self.single_atom_RO2_threshold = self.single_atom_threshold*self.t_SPCM_second_shot

        self.override_ExperimentVariables_dict = eval(self.override_ExperimentVariables)
        assert type(self.override_ExperimentVariables_dict) == dict, \
            "override_ExperimentVariables should be a python dictionary"

        for variable, value in self.override_ExperimentVariables_dict.items():
            assert hasattr(self, variable), (f"There is no ExperimentVariable " + variable +
                                             ". Did you mistype it?")

        self.variables_and_bounds = eval(self.variables_and_bounds)

        self.var_and_bounds_objects = []
        min_bounds = []
        max_bounds = []
        for var_and_bounds in self.variables_and_bounds:
            opt_var = OptimizerVariable(var_and_bounds, self)
            self.var_and_bounds_objects.append(opt_var)
            min_bounds.append(opt_var.min_bound)
            max_bounds.append(opt_var.max_bound)

        self.experiment_name = self.experiment_function
        self.experiment_function = lambda: eval(self.experiment_name)(self)

        self.cost_name = self.cost_function
        self.cost_function = lambda: eval(self.cost_name)(self)

        self.measurement = 0
        self.iteration = 0
        self.cost_dataset = "cost"
        self.initial_cost = 0
        self.best_cost = 0

        # instantiate the M-LOOP interface
        interface = MLOOPInterface()
        interface.get_next_cost_dict = self.get_next_cost_dict_for_mloop

        self.n_params = len(self.var_and_bounds_objects)

        self.mloop_controller = mlc.create_controller(interface,
                                                      max_num_runs=self.max_runs,
                                                      target_cost=self.target_cost,
                                                      num_params=self.n_params,
                                                      min_boundary=min_bounds,
                                                      max_boundary=max_bounds)

        self.measurement = 0
        self.SPCM0_RO1 = 0
        self.SPCM0_RO2 = 0

        # if there are multiple experiments in the schedule, then there might be something that has updated the datasets
        # e.g., as a result of an optimization scan. We want to make sure that this experiment uses the most up-to-date
        # datasets. However, ARTIQ runs build and prepare while the previous experiment is running, so our base.build
        # and base.prepare calls might use outdated datasets. This checks if there are other experiments scheduled
        # earlier than this one, and sets a flag so that base.build and base.prepare can be called again if needed.
        status_dict = self.scheduler.get_status()
        my_rid = self.scheduler.rid
        earlier_experiments = len([rid for rid, _ in status_dict.items() if rid < my_rid])
        logging.info("my rid is", my_rid, ", and there are", earlier_experiments,
                     " experiment(s) that I am waiting on to run")
        self.needs_fresh_build = earlier_experiments > 0

    def initialize_datasets(self):

        self.base.initialize_datasets()

        self.set_dataset(self.cost_dataset,
                         [0.0],
                         broadcast=True)

        self.set_dataset("optimizer_vars_dataset", [var.name for var in self.var_and_bounds_objects], broadcast=True)
        self.set_dataset("optimizer_bounds", [(var.min_bound, var.max_bound) for
                                                      var in self.var_and_bounds_objects], broadcast=True)
        self.optimizer_var_datasets = []

        for i in range(self.n_params):
            self.optimizer_var_datasets.append("optimizer_var"+str(i))
            self.set_dataset(self.optimizer_var_datasets[i], [0.0], broadcast=True)

        for var, val in self.override_ExperimentVariables_dict.items():
            self.set_dataset(var, val)

        value = 0.0
        for ch_i in range(len(self.laser_stabilizer.all_channels)):
            self.set_dataset(self.laser_stabilizer.all_channels[ch_i].dB_history_dataset,
                             [float(self.initial_RF_dB_values[ch_i])], broadcast=True, persist=True)


        # self.set_dataset('SPCM0_RO1_current_iteration', [0], broadcast=True)
        # self.set_dataset('SPCM0_RO2_current_iteration', [0], broadcast=True)
        # self.set_dataset('SPCM1_RO1_current_iteration', [0], broadcast=True)
        # self.set_dataset('SPCM1_RO2_current_iteration', [0], broadcast=True)


    def reset_datasets(self):
        """
        set datasets that are redefined each iteration.

        typically these datasets are used for plotting which would be meaningless if we continued to append to the SPCM0_RO1,
        e.g. for the second readout histogram which we expect in general will change as experiment parameters induce
        different amount of atom loss.
        :return:
        """
        self.set_dataset('SPCM0_RO1_current_iteration', [0], broadcast=True)
        self.set_dataset('SPCM0_RO2_current_iteration', [0], broadcast=True)
        self.set_dataset('SPCM1_RO1_current_iteration', [0], broadcast=True)
        self.set_dataset('SPCM1_RO2_current_iteration', [0], broadcast=True)
        self.set_dataset('BothSPCMs_RO1_current_iteration', [0], broadcast=True)
        self.set_dataset('BothSPCMs_RO2_current_iteration', [0], broadcast=True)

    def run(self):

        if self.needs_fresh_build:
            self.base.build()
            self.base.prepare()

        self.initialize_datasets()

        # override specific variables. this will apply to the entire scan, so it is outside the loops
        for variable, value in self.override_ExperimentVariables_dict.items():
            setattr(self, variable, value)

        # self.warm_up()
        self.initialize_hardware()  #todo: this is done in optimization routine

        # params is ignored when checking the initial cost
        cost = self.optimization_routine(params=[0.0]*self.n_params, check_initial_cost=True)
        self.mloop_controller.optimize()

        ### After optimiation ###
        print('Best parameters found:')
        print(self.mloop_controller.best_params)
        best_params = self.mloop_controller.best_params
        print(best_params)
        # self.set_experiment_variables_to_best_params(best_params)

        self.print_async("initial cost:", self.initial_cost)
        self.print_async("best cost:", self.best_cost)

        # todo: set best params to dataset

        # write the h5 file here in case worker refuses to die
        self.write_results({'name':self.experiment_name[:-11]+"_optimized_for_"+self.cost_name[:-5]})
        
    @kernel
    def initialize_hardware(self):
        self.base.initialize_hardware()

    @kernel
    def warm_up(self):
        """hardware init and turn things on"""

        self.core.reset()

        self.zotino0.set_dac([self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
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
        delay(2 * s)

        self.core.break_realtime()
        # warm up to get make sure we get to the setpoints
        for i in range(10):
            self.laser_stabilizer.run()
        self.dds_FORT.sw.on()

    def get_cost(self) -> TInt32:
        """cost function wrapper"""
        return self.cost_function()

    def optimization_routine(self, params: TArray(TFloat), check_initial_cost=False) -> TInt32:
        """
        the function that will be called by the optimizer.

        for use with M-LOOP, this should be called in the interface's "get_next_cost_dict"
        method.

        params: array of float values which we are trying to optimize
        check_initial_cost: boolean, default False. If true, evaluate the cost using the default values for the
            experiment variables we're trying to optimize. this will allow us to see if the optimizer has made
            any improvement.
        return:
            cost: the cost for the optimizer
        """

        # todo: set the experiment variables from params. see GeneralVariableScan

        for i in range(self.n_params):
            if not check_initial_cost:
                param_val = params[i]
                setattr(self, self.var_and_bounds_objects[i].name, param_val)
                self.append_to_dataset(self.optimizer_var_datasets[i], param_val)
            else:
                self.set_dataset(self.optimizer_var_datasets[i], [self.var_and_bounds_objects[i].default_value],
                                 broadcast=True)

        # todo: delete
        print("in opt. routine:", params)

        self.initialize_hardware()
        self.reset_datasets()

        # the measurement loop.
        self.experiment_function()

        cost = self.get_cost()
        # todo: delete
        # if check_initial_cost:
        #     print("this is the starting cost")
        # print("inside optimization routine:", self.iteration, cost)

        self.iteration += 1
        self.set_dataset("iteration", self.iteration, broadcast=True)

        if not check_initial_cost:
            self.append_to_dataset(self.cost_dataset, cost)
            if cost < self.best_cost:
                self.best_cost = cost
                self.print_async("new best cost:", self.best_cost)
                for i in range(self.n_params):
                    param_val = params[i]
                    self.print_async(self.var_and_bounds_objects[i].name, param_val)

                if self.set_best_parameters_at_finish:
                    for i in range(self.n_params):
                        param_val = params[i]
                        self.set_dataset(self.var_and_bounds_objects[i].name, float(param_val),
                                         broadcast=True,
                                         persist=True)
        else:
            self.set_dataset(self.cost_dataset, [cost], broadcast=True)
            self.initial_cost = cost

        # write the h5 file here in case we want to abort the experiment early
        self.write_results({'name': self.experiment_name[:-11] + "_optimized_for_" + self.cost_name[:-5]})

        return cost

    # this is fine for atom retention experiments, which is nearly
    # all experiments
    def get_next_cost_dict_for_mloop(self, params_dict):

        # Get parameters from the provided dictionary
        params = params_dict['params']

        cost = self.optimization_routine(params)
        # a proxy for the uncertainty, since the cost is typically -1*(number of atoms detected)
        uncertainty = 1/np.sqrt(-1*cost) if cost < 0 else 0

        cost_dict = {'cost': cost, 'uncer': uncertainty}
        return cost_dict

    @kernel  # shouldn't need to run this on the kernel
    def set_experiment_variables_to_best_params(self, best_params: TArray(TFloat)):
        self.core.reset()
        delay(1 * ms)

        if self.best_cost < self.initial_cost and self.set_best_parameters_at_finish:
            for i in range(self.n_params):
                param_val = best_params[i]
                self.set_dataset(self.var_and_bounds_objects[i].name, param_val,
                                 broadcast=True,
                                 persist=True)

            self.print_async("initial cost:", self.initial_cost)
            self.print_async("best cost:", self.best_cost)

    def analyze(self):
        mlv.show_all_default_visualizations(self.mloop_controller)


