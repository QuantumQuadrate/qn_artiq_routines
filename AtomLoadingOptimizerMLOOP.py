"""
Atom loading optimization which positions the MOT using M-LOOP
"""

from artiq.experiment import *
import numpy as np
import scipy as sp

#Imports for M-LOOP
import mloop.interfaces as mli
import mloop.controllers as mlc
import mloop.visualizations as mlv

from utilities.BaseExperiment import BaseExperiment


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
        # # Get parameters from the provided dictionary
        # params = params_dict['params']
        #
        # # Here you can include the code to run your experiment given a particular set of parameters
        # # In this example we will just evaluate a sum of sinc functions
        # cost = -np.sum(np.sinc(params - self.minimum_params))
        # # There is no uncertainty in our result
        # uncer = 0
        # # The evaluation will always be a success
        # bad = False
        # # Add a small time delay to mimic a real experiment
        # time.sleep(1)
        #
        # # The cost, uncertainty and bad boolean must all be returned as a dictionary
        # # You can include other variables you want to record as well if you want
        # cost_dict = {'cost': cost, 'uncer': uncer, 'bad': bad}
        # return cost_dict


class AtomLoadingOptimizerMLOOP(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        # overwrite the experiment variables of the same names
        # self.setattr_argument("t_SPCM_exposure", NumberValue(10 * ms, unit='ms'))
        self.setattr_argument("atom_counts_threshold", NumberValue(220))
        self.setattr_argument("t_MOT_loading", NumberValue(500 * ms, unit='ms'))
        self.setattr_argument("n_measurements", NumberValue(400, type='int', scale=1, ndecimals=0, step=1))
        self.setattr_argument("set_best_coil_volts_at_finish", BooleanValue(True))

        group = "differential coil volts boundaries"
        self.setattr_argument("use_differential_boundaries",BooleanValue(True),group)
        self.setattr_argument("dV_AZ_bottom", NumberValue(0.05, unit="V"), group)
        self.setattr_argument("dV_AZ_top", NumberValue(0.05, unit="V"), group)
        self.setattr_argument("dV_AX", NumberValue(0.05, unit="V"), group)
        self.setattr_argument("dV_AY", NumberValue(0.05, unit="V"), group)
        group = "absolute coil volts boundaries"
        self.setattr_argument("V_AZ_bottom_min", NumberValue(3.5,unit="V"), group)
        self.setattr_argument("V_AZ_bottom_max", NumberValue(4.8,unit="V"), group)
        self.setattr_argument("V_AZ_top_min", NumberValue(4.8,unit="V"), group)
        self.setattr_argument("V_AZ_top_max", NumberValue(4.8,unit="V"), group)
        self.setattr_argument("V_AX_min", NumberValue(-0.5,unit="V"), group)
        self.setattr_argument("V_AX_max", NumberValue(0.5,unit="V"), group)
        self.setattr_argument("V_AY_min", NumberValue(-0.5, unit="V"), group)
        self.setattr_argument("V_AY_max", NumberValue(0.5, unit="V"), group)

        group1 = "optimizer settings"
        self.setattr_argument("max_runs",NumberValue(70, type='int', scale=1, ndecimals=0, step=1),group1)

        # this should be close to the mean signal from the atom
        self.base.set_datasets_from_gui_args()
        print("build - done")

    def prepare(self):
        self.base.prepare()

        self.sampler_buffer = np.zeros(8)
        # self.minimize_args = eval(self.scipy_minimize_args)
        # self.minimize_options = eval(self.scipy_minimize_options)

        # whether to use boundaries which are defined as +/- dV around the current settings
        if self.use_differential_boundaries:
            # override the absolute boundaries
            self.V_AZ_bottom_min = self.AZ_bottom_volts_MOT - self.dV_AZ_bottom
            self.V_AZ_top_min = self.AZ_top_volts_MOT - self.dV_AZ_top
            self.V_AX_min = self.AX_volts_MOT - self.dV_AX
            self.V_AY_min = self.AY_volts_MOT - self.dV_AY

            self.V_AZ_bottom_max = self.AZ_bottom_volts_MOT + self.dV_AZ_bottom
            self.V_AZ_top_max = self.AZ_top_volts_MOT + self.dV_AZ_top
            self.V_AX_max = self.AX_volts_MOT + self.dV_AX
            self.V_AY_max = self.AY_volts_MOT + self.dV_AY

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

        # Next create the controller. Provide it with your interface and any options you want to set
        self.mloop_controller = mlc.create_controller(interface,
                                           max_num_runs=self.max_runs,
                                           target_cost=-0.5*self.n_measurements, # number of atoms to load
                                           num_params=4,
                                           min_boundary=[self.V_AZ_bottom_min,
                                                         self.V_AZ_top_min,
                                                         self.V_AX_min,
                                                         self.V_AY_min],
                                           max_boundary=[self.V_AZ_bottom_max,
                                                         self.V_AZ_top_max,
                                                         self.V_AX_max,
                                                         self.V_AY_max])

    def run(self):
        self.initialize_hardware()
        self.warm_up()

        self.mloop_controller.optimize()

        print('Best parameters found:')
        print(self.mloop_controller.best_params)
        best_volts = self.mloop_controller.best_params

        volt_datasets = ["AZ_bottom_volts_MOT", "AZ_top_volts_MOT", "AX_volts_MOT", "AY_volts_MOT"]
        if self.set_best_coil_volts_at_finish:
            print("updating coil values")
            for i in range(4):
                self.set_dataset(volt_datasets[i], float(best_volts[i]), broadcast=True, persist=True)

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
        self.laser_stabilizer.run()
        delay(1 * ms)
        self.dds_FORT.sw.on()
        delay(1*s)

    def get_cost(self, data: TArray(TFloat,1)) -> TInt32:
        q = 0 # whether the counts exceeded the atom threshold
        atoms_loaded = 0
        q_last = (data[0] > self.atom_counts_threshold)
        for x in data[1:]:
            q = x > self.atom_counts_threshold
            if q != q_last and q_last:
                atoms_loaded += 1
                # todo: replace with logging statement?
                # self.print_async("optimizer found the single atom signal!")
            q_last = q
        # return -1 if 0 < atoms_loaded < 10 else -1*atoms_loaded # threshold should depend on number of measurements
        return -1 * atoms_loaded

    @kernel
    def optimization_routine(self, coil_values: TArray(TFloat)) -> TInt32:
        """
        the function that will be called by the optimizer.

        for use with M-LOOP, this should be called in the interface's "get_next_cost_dict"
        method.

        coil_values: array of coil control volt values
        return:
            cost: the cost for the optimizer
        """
        self.core.reset()
        delay(1*ms)

        self.zotino0.set_dac(coil_values, channels=self.coil_channels)
        delay(1*ms)

        self.laser_stabilizer.run()
        delay(1*ms)
        self.dds_FORT.sw.on()

        delay(self.t_MOT_loading)

        # reset the counts dataset each run so we don't overwhelm the dashboard when plotting
        self.set_dataset(self.count_rate_dataset,[0.0],broadcast=True)

        for i in range(self.n_measurements):
            t_end = self.ttl0.gate_rising(self.t_SPCM_exposure)
            counts_per_s = self.ttl0.count(t_end)/self.t_SPCM_exposure
            delay(1 * ms)
            self.append_to_dataset(self.count_rate_dataset, counts_per_s)
            self.counts_list[i] = counts_per_s*self.t_SPCM_exposure

        cost = self.get_cost(self.counts_list)
        self.append_to_dataset(self.cost_dataset, cost)

        return cost

    def get_next_cost_dict_for_mloop(self,params_dict):

        # Get parameters from the provided dictionary
        params = params_dict['params']

        cost = self.optimization_routine(params)
        uncertainty = 1/np.sqrt(-1*cost) if cost < 0 else 0

        cost_dict = {'cost': cost, 'uncer': uncertainty}
        return cost_dict


