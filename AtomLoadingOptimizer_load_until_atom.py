"""
Atom loading optimization to optimize atom loading time with load_MOT_and_FORT_until_atom function.
MOT and FORT beams are turned on and the time to load an atom is minimized.

"""

from artiq.experiment import *
import numpy as np
import scipy as sp
from skimage.filters import threshold_otsu
import logging

### Imports for M-LOOP
import mloop.interfaces as mli
import mloop.controllers as mlc
import mloop.visualizations as mlv

from utilities.BaseExperiment import BaseExperiment


### Declare your custom class that inherits from the Interface class
class MLOOPInterface(mli.Interface):

    ### Initialization of the interface, including this method is optional
    def __init__(self):
        # You must include the super command to call the parent class, Interface, constructor
        super(MLOOPInterface, self).__init__()

    ### You must include the get_next_cost_dict method in your class
    ### this method is called whenever M-LOOP wants to run an experiment
    def get_next_cost_dict(self, params_dict):
        pass


class AtomLoadingOptimizer_load_until_atom(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        ### overwrite the experiment variables of the same names
        # self.setattr_argument("t_SPCM_exposure", NumberValue(10 * ms, unit='ms'))
        self.setattr_argument("atom_counts_per_s_threshold", NumberValue(14000))
        self.setattr_argument("n_measurements", NumberValue(10, type='int', scale=1, ndecimals=0, step=1))
        self.setattr_argument("set_best_parameters_at_finish", BooleanValue(True))
        self.both_mode = "coils and beam powers"
        self.beam_mode = "beam powers only"
        self.coil_mode = "coils only"
        self.setattr_argument("what_to_tune", EnumerationValue([self.both_mode, self.coil_mode, self.beam_mode]))

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

        group = "beam tuning settings"
        self.setattr_argument("max_set_point_percent_deviation_plus", NumberValue(0.1), group)
        self.setattr_argument("max_set_point_percent_deviation_minus", NumberValue(0.1), group)

        ### we can balance the z beams with confidence by measuring the powers outside the chamber,
        ### so unless we are fine tuning loading, we may want trust our initial manual balancing
        self.setattr_argument("disable_z_beam_tuning", BooleanValue(False), group)

        group1 = "optimizer settings"
        self.setattr_argument("max_runs",NumberValue(100, type='int', scale=1, ndecimals=0, step=1),group1)

        self.base.set_datasets_from_gui_args()
        print("build - done")

    def prepare(self):
        self.base.prepare()

        self.atom_counts_threshold = self.atom_counts_per_s_threshold*self.t_SPCM_exposure

        ### whether to use boundaries which are defined as +/- dV around the current settings
        if self.use_differential_boundaries:
            ### override the absolute boundaries
            self.V_AZ_bottom_min = self.AZ_bottom_volts_MOT - self.dV_AZ_bottom
            self.V_AZ_top_min = self.AZ_top_volts_MOT - self.dV_AZ_top
            self.V_AX_min = self.AX_volts_MOT - self.dV_AX
            self.V_AY_min = self.AY_volts_MOT - self.dV_AY

            self.V_AZ_bottom_max = self.AZ_bottom_volts_MOT + self.dV_AZ_bottom
            self.V_AZ_top_max = self.AZ_top_volts_MOT + self.dV_AZ_top
            self.V_AX_max = self.AX_volts_MOT + self.dV_AX
            self.V_AY_max = self.AY_volts_MOT + self.dV_AY

        self.coil_values = np.array([self.AZ_bottom_volts_MOT,
                                     self.AZ_top_volts_MOT,
                                     self.AX_volts_MOT,
                                     self.AY_volts_MOT])
        self.coil_values_RO = np.array([self.AZ_bottom_volts_PGC,
                                     - self.AZ_bottom_volts_PGC,
                                     self.AX_volts_PGC,
                                     self.AY_volts_PGC])

        self.volt_datasets = ["AZ_bottom_volts_MOT", "AZ_top_volts_MOT", "AX_volts_MOT", "AY_volts_MOT"]
        self.setpoint_datasets = ["set_point_PD1_AOM_A1", "set_point_PD2_AOM_A2", "set_point_PD3_AOM_A3",
                                  "set_point_PD4_AOM_A4", "set_point_PD5_AOM_A5", "set_point_PD6_AOM_A6"]
        self.default_setpoints = np.array([getattr(self, dataset) for dataset in self.setpoint_datasets])

        self.atom_loading_time_list = np.zeros(self.n_measurements)
        self.set_dataset(self.BothSPCMs_rate_dataset,
                         [0.0],
                         broadcast=True)

        self.cost_dataset = "cost"
        self.current_best_cost = 0
        self.set_dataset(self.cost_dataset,
                         [self.current_best_cost],
                         broadcast=True)

        ### instantiate the M-LOOP interface
        interface = MLOOPInterface()
        interface.get_next_cost_dict = self.get_next_cost_dict_for_mloop

        min_bounds = []
        max_bounds = []

        self.tune_beams = self.what_to_tune == self.beam_mode or self.what_to_tune == self.both_mode
        self.tune_coils = self.what_to_tune == self.coil_mode or self.what_to_tune == self.both_mode

        if self.tune_coils:
            print("MLOOP will tune coil volts")

            min_bounds += [self.V_AZ_bottom_min,
                           self.V_AZ_top_min,
                           self.V_AX_min,
                           self.V_AY_min]

            max_bounds += [self.V_AZ_bottom_max,
                           self.V_AZ_top_max,
                           self.V_AX_max,
                           self.V_AY_max]

        if self.tune_beams:
            print("MLOOP will tune beam powers")

            min_bounds += [1 - self.max_set_point_percent_deviation_minus,
                           1 - self.max_set_point_percent_deviation_minus,
                           1 - self.max_set_point_percent_deviation_minus,
                           1 - self.max_set_point_percent_deviation_minus]

            max_bounds += [1 + self.max_set_point_percent_deviation_plus,
                           1 + self.max_set_point_percent_deviation_plus,
                           1 + self.max_set_point_percent_deviation_plus,
                           1 + self.max_set_point_percent_deviation_plus]

            if not self.disable_z_beam_tuning:
                # append additional bounds for MOT5,6
                min_bounds.append(1 - self.max_set_point_percent_deviation_minus)
                min_bounds.append(1 - self.max_set_point_percent_deviation_minus)
                max_bounds.append(1 + self.max_set_point_percent_deviation_plus)
                max_bounds.append(1 + self.max_set_point_percent_deviation_plus)
        n_params = len(max_bounds)

        self.best_params = np.zeros(n_params)

        print("max bounds")
        print(max_bounds)
        print("min bounds")
        print(min_bounds)

        self.mloop_controller = mlc.create_controller(interface,
                                           max_num_runs=self.max_runs,
                                           target_cost=-3000.0, ### target atom_loading_time (in s) = -1000/target_cost.
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

    @kernel
    def warm_up(self):
        """hardware init and turn things on"""

        self.core.reset()

        delay(1*ms)

        ### Set the coils to MOT loading setting
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
            channels=self.coil_channels)

        ### set the cooling DP AOM to the MOT settings
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)

        self.dds_cooling_DP.sw.on()  ### turn on cooling
        self.ttl_repump_switch.off()  ### turn on MOT RP

        self.dds_AOM_A1.sw.on()
        self.dds_AOM_A2.sw.on()
        self.dds_AOM_A3.sw.on()
        self.dds_AOM_A4.sw.on()
        delay(0.1 * ms)
        self.dds_AOM_A5.sw.on()
        self.dds_AOM_A6.sw.on()

        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude)
        self.dds_FORT.sw.on()

        ### delay for AOMs to thermalize
        delay(500 * ms)

        # self.core.break_realtime()
        ### warm up to make sure we get to the setpoints
        for i in range(10):
            self.laser_stabilizer.run()

        ### Turning off AOMs to be ready to start atom loading from scratch
        self.dds_FORT.sw.off() ### turn off FORT
        self.dds_cooling_DP.sw.off()  ### turn off cooling
        self.ttl_repump_switch.on()  ### turn off MOT RP


    def get_cost(self, data: TArray(TFloat,1)) -> TInt32:
        total_t = 0.0
        for t in data:
            total_t += -1000.0 / t
        average_t = total_t / len(data) ### Though I am naming these at _t, these are indeed 1/t to calculate the cost
        return int(round(average_t))


    @kernel
    def optimization_routine(self, params: TArray(TFloat)) -> TInt32:
        """
        For use with M-LOOP, this should be called in the interface's "get_next_cost_dict"
        method.

        params: array of float values which we are trying to optimize
        return:
            cost: the cost for the optimizer
        """

        self.core.reset()
        delay(1*ms)

        # self.zotino0.set_dac([3.5], self.Osc_trig_channel)  ### for triggering oscilloscope
        # delay(0.1 * ms)
        # self.zotino0.set_dac([0.0], self.Osc_trig_channel)

        if self.tune_coils:
            self.coil_values = params[:4]
            if self.tune_beams:
                setpoint_multipliers = params[4:]
            else:
                setpoint_multipliers = np.array([1.0,1.0,1.0,1.0,1.0,1.0])
        else:
            setpoint_multipliers = params

        if self.tune_beams:
            self.stabilizer_AOM_A1.set_points[0] = self.default_setpoints[0] * setpoint_multipliers[0]
            self.stabilizer_AOM_A2.set_points[0] = self.default_setpoints[1] * setpoint_multipliers[1]
            self.stabilizer_AOM_A3.set_points[0] = self.default_setpoints[2] * setpoint_multipliers[2]
            self.stabilizer_AOM_A4.set_points[0] = self.default_setpoints[3] * setpoint_multipliers[3]
            if not self.disable_z_beam_tuning and not self.what_to_tune == self.coil_mode:
                self.stabilizer_AOM_A5.set_points[0] = self.default_setpoints[4] * setpoint_multipliers[4]
                self.stabilizer_AOM_A6.set_points[0] = self.default_setpoints[5] * setpoint_multipliers[5]

            for i in range(3):
                self.laser_stabilizer.run()
        else:
            self.laser_stabilizer.run()

        if self.tune_coils:
            self.zotino0.set_dac(self.coil_values, channels=self.coil_channels)

        delay(1 * ms)

        ##################### This is the core of the optimizer that runs the sequence and get a cost:

        ### reset the counts dataset each run so we don't overwhelm the dashboard when plotting
        self.set_dataset(self.BothSPCMs_rate_dataset, [0.0], broadcast=True)

        self.laser_stabilizer.run()

        for i in range(self.n_measurements):
            delay(1 * ms)
            self.dds_cooling_DP.sw.on()  ### turn on cooling
            self.ttl_repump_switch.off()  ### turn on MOT RP

            self.dds_AOM_A1.sw.on()
            self.dds_AOM_A2.sw.on()
            self.dds_AOM_A3.sw.on()
            self.dds_AOM_A4.sw.on()
            delay(0.1 * ms)
            self.dds_AOM_A5.sw.on()
            self.dds_AOM_A6.sw.on()
            self.dds_FORT.sw.on()

            delay(1 * ms)
            # self.zotino0.set_dac([3.5], self.UV_trig_channel) ### for some reason it complains about this line. So, no UV for now
            ### in the optimizer.
            delay(1*ms)
            # self.zotino0.set_dac([3.5], self.Osc_trig_channel)  ### for triggering oscilloscope

            max_tries = 100  ### Maximum number of attempts before running the feedback
            atom_check_time   = 20 * ms
            atom_loaded = False
            try_n = 0
            t_before_atom = now_mu()  ### is used to calculate the loading time of atoms by atom_loading_time = t_after_atom - t_before_atom
            t_after_atom = now_mu()

            while not atom_loaded and try_n < max_tries:
                delay(100 * us)  ### Needs a delay of about 100us or maybe less
                if self.which_node == 'alice':
                    with parallel:
                        self.ttl_SPCM0_counter.gate_rising(atom_check_time)
                        self.ttl_SPCM1_counter.gate_rising(atom_check_time)

                    BothSPCMs_atom_check = int((self.ttl_SPCM0_counter.fetch_count() + self.ttl_SPCM1_counter.fetch_count()) / 2)
                else:
                    with parallel:
                        self.ttl_SPCM0_counter.gate_rising(atom_check_time)

                    BothSPCMs_atom_check = int(self.ttl_SPCM0_counter.fetch_count())

                BothSPCMs_counts_per_s = BothSPCMs_atom_check / atom_check_time
                delay(1 * ms)
                self.append_to_dataset(self.BothSPCMs_rate_dataset, BothSPCMs_counts_per_s)
                try_n += 1

                if BothSPCMs_counts_per_s > self.atom_counts_per_s_threshold:
                    delay(100 * us)  ### Needs a delay of about 100us or maybe less
                    atom_loaded = True

            if atom_loaded:
                t_after_atom = now_mu()
                atom_loading_time = self.core.mu_to_seconds(t_after_atom - t_before_atom)
            else:
                atom_loading_time = 10e9 ### Just a large number to show no atom loading. This is compared to the typical 0.5 second atom loading.

            # self.zotino0.set_dac([0.0], self.UV_trig_channel)
            delay(100*us)
            # self.zotino0.set_dac([0.0], self.Osc_trig_channel)
            self.atom_loading_time_list[i] = atom_loading_time


            delay(1 * ms)
            ### Turning off AOMs to be ready to start atom loading from scratch
            self.ttl_repump_switch.on()  ### turn off MOT RP
            self.dds_cooling_DP.sw.off()  ### turn off cooling
            self.dds_FORT.sw.off()  ### turn off FORT
            self.dds_AOM_A1.sw.on()
            self.dds_AOM_A2.sw.on()
            delay(0.1 * ms)
            self.dds_AOM_A3.sw.on()
            self.dds_AOM_A4.sw.on()
            self.dds_AOM_A5.sw.on()
            self.dds_AOM_A6.sw.on()
            delay(300 * ms)  ### to dissipate MOT


        cost = self.get_cost(self.atom_loading_time_list)
        self.append_to_dataset(self.cost_dataset, cost)

        ################################################################################

        param_idx = 0
        if cost < self.current_best_cost:
            self.current_best_cost = cost
            self.print_async("NEW BEST COST:", cost)
            if self.tune_coils:
                self.print_async("BEST coil values:", params[:4])
                self.best_params[:4] = params[:4]
                param_idx = 3
            if self.tune_beams:
                for i in range(4):
                    self.best_params[param_idx + i] = self.default_setpoints[i] * setpoint_multipliers[i]
                    self.print_async("BEST setpoint",i+1,self.best_params[param_idx + i])
                if not self.disable_z_beam_tuning:
                    self.best_params[param_idx + 4] = self.default_setpoints[4] * setpoint_multipliers[4]
                    self.best_params[param_idx + 5] = self.default_setpoints[5] * setpoint_multipliers[5]
                    self.print_async("BEST setpoint", 5, self.best_params[param_idx + 4])
                    self.print_async("BEST setpoint", 6, self.best_params[param_idx + 5])
            self.set_dataset("best params", self.best_params)

        return cost



    @kernel
    def optimization_routine_test(self, params: TArray(TFloat)) -> TInt32:
        """
        Added by Eunji.
        Requires explanation.
        Does not turn on the coils in mode "beam powers only".


        For use with M-LOOP, this should be called in the interface's "get_next_cost_dict"
        method.

        Changes:
        * every measurement starts by setting MOT_coils.
        * then it ends wit RO_coil settings.



        params: array of float values which we are trying to optimize
        return:
            cost: the cost for the optimizer
        """

        self.core.reset()
        delay(1*ms)

        # self.zotino0.set_dac([3.5], self.Osc_trig_channel)  ### for triggering oscilloscope
        # delay(0.1 * ms)
        # self.zotino0.set_dac([0.0], self.Osc_trig_channel)

        if self.tune_coils:
            self.coil_values = params[:4]
            if self.tune_beams:
                setpoint_multipliers = params[4:]
            else:
                setpoint_multipliers = np.array([1.0,1.0,1.0,1.0,1.0,1.0])
        else:
            setpoint_multipliers = params

        if self.tune_beams:
            self.stabilizer_AOM_A1.set_points[0] = self.default_setpoints[0] * setpoint_multipliers[0]
            self.stabilizer_AOM_A2.set_points[0] = self.default_setpoints[1] * setpoint_multipliers[1]
            self.stabilizer_AOM_A3.set_points[0] = self.default_setpoints[2] * setpoint_multipliers[2]
            self.stabilizer_AOM_A4.set_points[0] = self.default_setpoints[3] * setpoint_multipliers[3]
            if not self.disable_z_beam_tuning and not self.what_to_tune == self.coil_mode:
                self.stabilizer_AOM_A5.set_points[0] = self.default_setpoints[4] * setpoint_multipliers[4]
                self.stabilizer_AOM_A6.set_points[0] = self.default_setpoints[5] * setpoint_multipliers[5]

            for i in range(3):
                self.laser_stabilizer.run()
        else:
            self.laser_stabilizer.run()

        # if self.tune_coils:
        #     self.zotino0.set_dac(self.coil_values, channels=self.coil_channels)

        delay(1 * ms)

        ##################### This is the core of the optimizer that runs the sequence and get a cost:

        ### reset the counts dataset each run so we don't overwhelm the dashboard when plotting
        self.set_dataset(self.BothSPCMs_rate_dataset, [0.0], broadcast=True)

        for i in range(self.n_measurements):
            if self.tune_coils:
                self.zotino0.set_dac(self.coil_values, channels=self.coil_channels)
            delay(10*ms)
            # delay(1 * ms)
            self.laser_stabilizer.run()
            self.dds_cooling_DP.sw.on()  ### turn on cooling
            self.ttl_repump_switch.off()  ### turn on MOT RP

            self.dds_AOM_A1.sw.on()
            self.dds_AOM_A2.sw.on()
            self.dds_AOM_A3.sw.on()
            self.dds_AOM_A4.sw.on()
            delay(0.1 * ms)
            self.dds_AOM_A5.sw.on()
            self.dds_AOM_A6.sw.on()
            self.dds_FORT.sw.on()

            delay(1 * ms)

            max_tries = 100  ### Maximum number of attempts before running the feedback
            atom_check_time   = 20 * ms
            atom_loaded = False
            try_n = 0
            t_before_atom = now_mu()  ### is used to calculate the loading time of atoms by atom_loading_time = t_after_atom - t_before_atom
            t_after_atom = now_mu()

            while not atom_loaded and try_n < max_tries:
                delay(100 * us)  ### Needs a delay of about 100us or maybe less
                if self.which_node == 'alice':
                    with parallel:
                        self.ttl_SPCM0_counter.gate_rising(atom_check_time)
                        self.ttl_SPCM1_counter.gate_rising(atom_check_time)

                    BothSPCMs_atom_check = int((self.ttl_SPCM0_counter.fetch_count() + self.ttl_SPCM1_counter.fetch_count()) / 2)
                else:
                    with parallel:
                        self.ttl_SPCM0_counter.gate_rising(atom_check_time)

                    BothSPCMs_atom_check = int(self.ttl_SPCM0_counter.fetch_count())

                BothSPCMs_counts_per_s = BothSPCMs_atom_check / atom_check_time
                delay(1 * ms)
                self.append_to_dataset(self.BothSPCMs_rate_dataset, BothSPCMs_counts_per_s)
                try_n += 1

                if BothSPCMs_counts_per_s > self.atom_counts_per_s_threshold:
                    delay(100 * us)  ### Needs a delay of about 100us or maybe less
                    atom_loaded = True

            if atom_loaded:
                t_after_atom = now_mu()
                atom_loading_time = self.core.mu_to_seconds(t_after_atom - t_before_atom)
            else:
                atom_loading_time = 10e9 ### Just a large number to show no atom loading. This is compared to the typical 0.5 second atom loading.

            self.atom_loading_time_list[i] = atom_loading_time

            ##########checking
            self.zotino0.set_dac(self.coil_values_RO, channels=self.coil_channels)


            delay(1 * ms)
            ### Turning off AOMs to be ready to start atom loading from scratch
            self.ttl_repump_switch.on()  ### turn off MOT RP
            self.dds_cooling_DP.sw.off()  ### turn off cooling
            self.dds_FORT.sw.off()  ### turn off FORT
            delay(100 * ms)  ### to dissipate MOT

        # #############todo: why does it give me error?
        # # # delay(1 * ms)
        # # # ### set the coils to the readout settings
        # # # self.zotino0.set_dac(
        # # #     [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
        # # #     channels=self.coil_channels)
        # # delay(1 * ms)
        # # delay(1*s)
        # # if True:
        # #     ### set the coils to the readout settings
        # #     self.zotino0.set_dac([0.0*V,0.0*V,0.0*V,0.0*V],channels=self.coil_channels)
        #
        # self.zotino0.set_dac(self.coil_values_RO, channels=self.coil_channels)

        cost = self.get_cost(self.atom_loading_time_list)
        self.append_to_dataset(self.cost_dataset, cost)

        ################################################################################

        param_idx = 0
        if cost < self.current_best_cost:
            self.current_best_cost = cost
            self.print_async("NEW BEST COST:", cost)
            if self.tune_coils:
                self.print_async("BEST coil values:", params[:4])
                self.best_params[:4] = params[:4]
                param_idx = 3
            if self.tune_beams:
                for i in range(4):
                    self.best_params[param_idx + i] = self.default_setpoints[i] * setpoint_multipliers[i]
                    self.print_async("BEST setpoint",i+1,self.best_params[param_idx + i])
                if not self.disable_z_beam_tuning:
                    self.best_params[param_idx + 4] = self.default_setpoints[4] * setpoint_multipliers[4]
                    self.best_params[param_idx + 5] = self.default_setpoints[5] * setpoint_multipliers[5]
                    self.print_async("BEST setpoint", 5, self.best_params[param_idx + 4])
                    self.print_async("BEST setpoint", 6, self.best_params[param_idx + 5])
            self.set_dataset("best params", self.best_params)

        return cost


    def get_next_cost_dict_for_mloop(self,params_dict):

        ### Get parameters from the provided dictionary
        params = params_dict['params']

        cost = self.optimization_routine(params)
        # cost = self.optimization_routine_test(params) ### does not turn on coils in "beam powers only" mode.

        uncertainty = 1/np.sqrt(-1*cost) if cost < 0 else 0

        cost_dict = {'cost': cost, 'uncer': uncertainty}
        return cost_dict

    @kernel
    def set_experiment_variables_to_best_params(self, best_params: TArray(TFloat)):
        self.core.reset()
        delay(1 * ms)

        best_volts = self.coil_values
        best_setpoint_multipliers = self.default_setpoints

        if self.tune_coils:
            best_volts = best_params[:4]
            if self.tune_beams:
                best_setpoint_multipliers = best_params[4:]
        else:
            best_setpoint_multipliers = best_params

        if self.set_best_parameters_at_finish:
            if self.tune_coils:
                self.print_async("updating coil values")
                for i in range(4):
                    self.set_dataset(self.volt_datasets[i], float(best_volts[i]), broadcast=True, persist=True)
            if self.tune_beams:
                self.print_async("updating MOT beam setpoints")
                if self.disable_z_beam_tuning:
                    n_beams = 4
                else:
                    n_beams = 6

                for i in range(n_beams):
                    self.set_dataset(self.setpoint_datasets[i], self.default_setpoints[i] * best_setpoint_multipliers[i],
                                     broadcast=True,
                                     persist=True)

    def analyze(self):
        mlv.show_all_default_visualizations(self.mloop_controller)


