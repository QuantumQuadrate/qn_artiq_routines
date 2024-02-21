"""
For optimizing an experiment function over any number of defined variables

Any experiment function defined in experiment_functions.py can be optimized
over any set of ExperimentVariables

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
import subroutines.experiment_functions as exp_functions
import subroutines.cost_functions as cost_functions

class OptimizerVariable:
    """
    a container for the variables that we want to optimize.
    for streamlining setup
    """
    def __init__(self, var_and_bounds_tuple, experiment):
        self.var_name = var_and_bounds_tuple[0]
        self.is_absolute = 0 if var_and_bounds_tuple[3] == 'diff' else 1
        self.default_value = getattr(self,self.var_name)
        self.min_bound = var_and_bounds_tuple[1] + self.is_absolute*self.default_value
        self.max_bound = var_and_bounds_tuple[2] + self.is_absolute*self.default_value

# Declare your custom class that inherits from the Interface class
class MLOOPInterface(mli.Interface):

    # Initialization of the interface, including this method is optional
    def __init__(self):
        # You must include the super command to call the parent class, Interface, constructor
        super(MLOOPInterface, self).__init__()

        # Attributes of the interface can be added here
        # If you want to precalculate any variables etc. this is the place to do it
        # In this example we will just define the location of the minimum
        self.minimum_params = np.array([0, 0.1, -0.1])

    # You must include the get_next_cost_dict method in your class
    # this method is called whenever M-LOOP wants to run an experiment
    def get_next_cost_dict(self, params_dict):
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
            "[('dummy_variable',0*ms,100*ms,'abs'),(dummy_variable,-0.5*V,0.5*V,'diff')]"))

        # todo: think carefully about how to handle feeding back to things
        #  that are not parents of experiment directly, e.g. the setpoints
        #  for the MOT beams. I had to explicitly code this in in AtomLoadingOptimizer

        experiment_function_names_list = [x for x in dir(exp_functions)
                                          if ('__' not in x and str(
                type(getattr(exp_functions, x))) == "<class 'function'>"
                                              and 'experiment' in x)]

        cost_function_names_list = [x for x in dir(cost_functions)
                                          if ('__' not in x and str(
                type(getattr(cost_functions, x))) == "<class 'function'>"
                                              and 'cost' in x)]

        group1 = "optimizer settings"
        self.setattr_argument("max_runs",NumberValue(70, type='int', scale=1, ndecimals=0, step=1),group1)

        # todo: add keyword arguments for setting up M-LOOP
        #  there is a lot more available then what we've been using

        self.base.set_datasets_from_gui_args()
        logging.debug("build - done")

    def prepare(self):
        self.base.prepare()

        # todo: see GeneralVariableScan
        #  define the experiment and cost function.

        self.variables_and_bounds = eval(self.variables_and_bounds)
        var_and_bounds_objects = []
        min_bounds = []
        max_bounds = []
        for var_and_bounds in self.variables_and_bounds:
            opt_var = OptimizerVariable(var_and_bounds,self)
            var_and_bounds_objects.append(opt_var)
            min_bounds.append(var_and_bounds.min_bound)
            max_bounds.append(var_and_bounds.max_bound)

        # todo
        # 1. build the min and max boundaries list from the variables and bounds list
        #    this has to be built depending on whether the bounds are diff


        self.default_variable_values = [getattr(self, dataset) for dataset in self.setpoint_datasets]

        self.counts_list = np.zeros(self.n_measurements)
        self.set_dataset(self.count_rate_dataset,
                         [0.0],
                         broadcast=True)

        self.cost_dataset = "cost"
        self.set_dataset(self.cost_dataset,
                         [0.0],
                         broadcast=True)

        # instantiate the MLOOP interface
        interface = MLOOPInterface()
        interface.get_next_cost_dict = self.get_next_cost_dict_for_mloop

        n_params = len(var_and_bounds_objects)

        self.mloop_controller = mlc.create_controller(interface,
                                           max_num_runs=self.max_runs,
                                           target_cost=-0.5*self.n_measurements, # -1 * number of atoms to load
                                           num_params=n_params,
                                           min_boundary=min_bounds,
                                           max_boundary=max_bounds)

    def run(self):
        self.initialize_hardware()
        self.warm_up()

        self.mloop_controller.optimize()

        print('Best parameters found:')
        print(self.mloop_controller.best_params)
        best_params = self.mloop_controller.best_params
        self.set_experiment_variables_to_best_params(best_params)

    @kernel
    def initialize_hardware(self):
        self.base.initialize_hardware()

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

        # call the experiment function here


        cost = self.get_cost(self.counts_list)
        self.append_to_dataset(self.cost_dataset, cost)

        return cost

    # this is fine for atom retention experiments, which is nearly
    # all experiments
    def get_next_cost_dict_for_mloop(self,params_dict):

        # Get parameters from the provided dictionary
        params = params_dict['params']

        cost = self.optimization_routine(params)
        uncertainty = 1/np.sqrt(-1*cost) if cost < 0 else 0

        cost_dict = {'cost': cost, 'uncer': uncertainty}
        return cost_dict

    @kernel # shouldn't need to run this on the kernel
    def set_experiment_variables_to_best_params(self, best_params: TArray(TFloat)):
        # self.core.reset()
        # delay(1 * ms)

        # todo: generalize this.

        for i in range(n_beams):
            self.set_dataset(self.setpoint_datasets[i], self.default_setpoints[i] * best_setpoint_multipliers[i],
                             broadcast=True,
                             persist=True)

    def analyze(self):
        mlv.show_all_default_visualizations(self.mloop_controller)


