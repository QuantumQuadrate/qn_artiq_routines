"""
For optimizing an experiment function over any number of defined variables

Any experiment function defined in experiment_functions.py can be optimized
over any set of ExperimentVariables

Example strings for the variables_and_bounds argument:
    [('AZ_bottom_volts_RO',-0.4*V,0.0*V,'diff'),('AZ_top_volts_RO',-0.1*V,0.3*V,'diff'),
('AX_volts_RO',-0.1*V,0.3*V,'diff'),('AY_volts_RO',-0.1*V,0.3*V,'diff')]
    [('f_cooling_DP_RO',-2*MHz,2*MHz,'diff'),('p_cooling_DP_RO',0.75,1.0,'abs')]
    [('p_cooling_DP_PGC',0.1,1.0,'abs'),('f_cooling_DP_PGC',115*MHz,125*MHz,'abs'),
('t_PGC_in_MOT',1*ms,40*ms,'abs')]

    [('set_point_PD1_AOM_A1',0.5,1.1,'perc'),('set_point_PD2_AOM_A2',0.5,1.1,'perc'),
('set_point_PD3_AOM_A3',0.5,1.1,'perc'),('set_point_PD4_AOM_A4',0.5,1.1,'perc'),
('set_point_PD5_AOM_A5',0.5,1.1,'perc'),('set_point_PD6_AOM_A6',0.5,1.1,'perc'),
('AZ_bottom_MOT',-0.2*V,0.2*V,'diff'),('AZ_top_MOT',-0.2*V,0.2*V,'diff'),
('AX_volts_MOT',-0.2*V,0.2*V,'diff'),('AY_volts_MOT',-0.2*V,0.2*V,'diff')]

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
from subroutines.aom_feedback import AOMPowerStabilizer


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
    def get_next_cost_dict(self, params_dict): # we'll redefine this function later
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

        group1 = "optimizer settings"
        self.setattr_argument("max_runs", NumberValue(70, type='int', scale=1, ndecimals=0, step=1), group1)
        self.setattr_argument("target_cost", NumberValue(-100, type='int', scale=1, ndecimals=0, step=1), group1)

        # todo: add keyword arguments for setting up M-LOOP
        #  there is a lot more available than what we've been using

        self.base.set_datasets_from_gui_args()
        logging.debug("build - done")

    def prepare(self):
        self.base.prepare()

        self.atom_counts_threshold = self.single_atom_counts_per_s*self.t_SPCM_first_shot
        self.atom_counts2_threshold = self.single_atom_counts_per_s*self.t_SPCM_second_shot

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
        self.set_dataset(self.cost_dataset,
                         [0.0],
                         broadcast=True)

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

    def initialize_datasets(self):
        self.set_dataset("n_measurements", self.n_measurements, broadcast=True)
        self.set_dataset("iteration", 0, broadcast=True)

        self.set_dataset("photocounts", [0], broadcast=True)
        self.set_dataset("photocounts2", [0], broadcast=True)
        self.set_dataset("photocount_bins", [50], broadcast=True)

        self.set_dataset("optimizer_vars_dataset", [var.name for var in self.var_and_bounds_objects], broadcast=True)
        self.set_dataset("optimizer_bounds_dataset", [(var.min_bound, var.max_bound) for
                                                      var in self.var_and_bounds_objects], broadcast=True)
        self.optimizer_var_datasets = []

        for i in range(self.n_params):
            self.optimizer_var_datasets.append("optimizer_var"+str(i))
            self.set_dataset(self.optimizer_var_datasets[i], [0.0], broadcast=True)

    def reset_datasets(self):
        """
        set datasets that are redefined each iteration.

        typically these datasets are used for plotting which would be meaningless if we continued to append to the photocounts,
        e.g. for the second readout histogram which we expect in general will change as experiment parameters induce
        different amount of atom loss.
        :return:
        """
        self.set_dataset('photocounts_current_iteration', [0], broadcast=True)
        self.set_dataset('photocounts2_current_iteration', [0], broadcast=True)

    def run(self):
        self.initialize_datasets()
        self.initialize_hardware()
        self.warm_up()

        self.mloop_controller.optimize()

        print('Best parameters found:')
        print(self.mloop_controller.best_params)
        best_params = self.mloop_controller.best_params
        self.set_experiment_variables_to_best_params(best_params)

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

    # todo should check if the cost has to be int
    def optimization_routine(self, params: TArray(TFloat)) -> TInt32:
        """
        the function that will be called by the optimizer.

        for use with M-LOOP, this should be called in the interface's "get_next_cost_dict"
        method.

        params: array of float values which we are trying to optimize
        return:
            cost: the cost for the optimizer
        """

        # todo: set the experiment variables from params. see GeneralVariableScan

        for i in range(self.n_params):
            param_val = params[i]
            setattr(self, self.var_and_bounds_objects[i].name, param_val)
            # todo broadcast the values so we can see what's going on
            self.append_to_dataset(self.optimizer_var_datasets[i], param_val)

        self.initialize_hardware()
        self.reset_datasets()

        # the measurement loop.
        self.experiment_function()

        self.iteration += 1
        self.set_dataset("iteration", self.iteration, broadcast=True)

        cost = self.get_cost()
        self.append_to_dataset(self.cost_dataset, cost)

        return cost

    # this is fine for atom retention experiments, which is nearly
    # all experiments
    def get_next_cost_dict_for_mloop(self,params_dict):

        # Get parameters from the provided dictionary
        params = params_dict['params']

        cost = self.optimization_routine(params)
        # a proxy for the uncertainty, since the cost is typically -1*(number of atoms detected)
        uncertainty = 1/np.sqrt(-1*cost) if cost < 0 else 0

        cost_dict = {'cost': cost, 'uncer': uncertainty}
        return cost_dict

    @kernel # shouldn't need to run this on the kernel
    def set_experiment_variables_to_best_params(self, best_params: TArray(TFloat)):
        self.core.reset()
        delay(1 * ms)

        for i in range(self.n_params):
            param_val = best_params[i]
            self.set_dataset(self.var_and_bounds_objects[i].name, param_val,
                             broadcast=True,
                             persist=True)

    def analyze(self):
        mlv.show_all_default_visualizations(self.mloop_controller)


