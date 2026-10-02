"""Atom loading optimization for the master-satellite architecture.

Ported from standalone/AtomLoadingOptimizer_load_until_atom.py, which stays
untouched and keeps serving the standalone systems. MOT and FORT beams are
turned on and M-LOOP minimizes the time taken to load an atom, tuning the MOT
coil volts and/or the six MOT beam set points.

which_node selects Node1 or Node2 and is the only added GUI argument; two_nodes
is deliberately not supported, because atom loading tunes one node's coils and
beams against that node's own SPCM rate, and master-satellite laser feedback is
single-node only for now (BaseExperiment_master_satellite.prepare_laser_
stabilizer raises otherwise).

WHY THE SEQUENCE ITSELF NEEDED NO CHANGES

The four SPCMs are master-local on both nodes and Base publishes the legacy
aliases this file already used (LEGACY_SPCM_ALIASES):

    ttl_SPCM0 -> SPCM_H1        ttl_SPCM0_OtherNode -> SPCM_H2
    ttl_SPCM1 -> SPCM_V1        ttl_SPCM1_OtherNode -> SPCM_V2

and AllSPCMs_rate_dataset alongside them, so the counting loop, the cost
function and both optimization routines are byte-for-byte the standalone ones.
Coils, DDS aliases, the zotino and the stabilizers all arrive under their
ordinary legacy names through the selected node's device-alias projection.

WHAT DID CHANGE, AND WHY

* Base is BaseExperimentMasterSatellite, configured in single_node mode, and
  the class mixes in _DatasetRedirectMixin so hard-coded legacy dataset names
  resolve to the selected node.

* The best-parameter writes are resolved explicitly in prepare(). This is the
  one place a verbatim copy would have gone quietly wrong. Both lists are
  written with persist=True at the end of a successful optimization, but only
  the coil volts are covered by the dataset redirect (they are exactly
  COIL_CALIBRATION_DATASETS). The six set_point_PD*_AOM_A* names are NOT: the
  feedback map holds each channel's monitor and power datasets
  (MOT1_monitor, p_AOM_A1, ...), not its set point. Written unsuffixed they
  would land on names no master-satellite code reads, so the beam-power half
  of every optimization would be silently discarded while
  set_point_PD1_AOM_A1_<node> kept its old value. Both lists therefore go
  through resolve_experiment_variable_target, which is the resolver for
  experiment VARIABLES, and the redirect passes the already-suffixed names
  through unchanged.

* base.set_datasets_from_gui_args() has no master-satellite equivalent and is
  dropped. In the standalone stack it is what makes a GUI argument win over
  the persistent variable of the same name: it writes the submitted value to
  the dataset in build(), so prepare()'s dataset load reads it back. Here the
  two arguments that shadow persistent variables -- t_SPCM_exposure (node
  suffixed) and n_measurements (global) -- are captured before
  configure_execution and re-asserted after, the same way MonitorSPCMinApplet
  and the GVS mixin do it. ARTIQ archives every submitted argument in
  expid["arguments"] regardless, which is where the analysis notebooks read
  them from.

* run() forces the idle node's beams, coils and microwave/RF DDS off before
  warming up, as every other single-node master-satellite experiment does.

* Every laser_stabilizer.run() is FOLLOWED by core.break_realtime(). This is
  the one place the sequence had to change, and it is not cosmetic: the
  master-satellite feedback variables are far heavier than the standalone
  ones -- aom_feedback_averages 20 vs 10, and aom_feedback_iterations capped
  at 200 (Node1) / 100 (Node2) vs 15 -- while subroutines/aom_feedback.py
  advances the timeline by only 0.2 ms per iteration and manages no slack of
  its own. Measured on hardware, one run returns roughly 5.8 ms behind the
  wall clock, so the next RTIO event underflows regardless of what it is:
  first attempt failed inside the warm-up loop at -4.9 ms, and re-arming
  BEFORE each run merely moved it to dds_FORT.sw.off() after the loop at
  -5.8 ms. After each run is the placement that covers both. Every one of
  these calls sits outside any measured interval (warm-up measures nothing;
  the others all precede t_before_atom), so moving the cursor changes no
  timing that matters.

* The transient plot datasets "cost" and "best params" stay unsuffixed, as in
  the standalone file. Only one node runs at a time, so there is nothing to
  collide with, and the applets read these names as they are.
"""

from artiq.experiment import *
import numpy as np
import logging

### Imports for M-LOOP
import mloop.interfaces as mli
import mloop.controllers as mlc
import mloop.visualizations as mlv

### Imported by name rather than with the `import *` MicrowaveScanOptimizer
### uses: this file needs exactly one function from that 9,500-line module, and
### naming it keeps the dependency visible.
from subroutines.experiment_functions import (
    run_feedback_and_record_FORT_MM_power,
)
from utilities.BaseExperiment_master_satellite import (
    BaseExperimentMasterSatellite,
    _DatasetRedirectMixin,
)


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


class AtomLoadingOptimizer_load_until_atom_master_satellite(
    _DatasetRedirectMixin, EnvExperiment
):
    """AtomLoadingOptimizer_load_until_atom_master_satellite

    Minimize the selected node's atom loading time with M-LOOP, tuning its MOT
    coil volts and/or MOT beam set points.
    """

    VALID_NODES = ("Node1", "Node2")
    # Tests inject a fake stabilizer class here; None selects the real
    # AOMPowerStabilizer from subroutines/aom_feedback.py.
    _stabilizer_factory = None

    ### The authoritative variables this experiment can persist at finish.
    ### Resolved to their node-suffixed names in prepare(); see the module
    ### docstring for why the set points cannot rely on the dataset redirect.
    COIL_VOLT_VARIABLES = (
        "AZ_bottom_volts_MOT",
        "AZ_top_volts_MOT",
        "AX_volts_MOT",
        "AY_volts_MOT",
    )
    SET_POINT_VARIABLES = (
        "set_point_PD1_AOM_A1",
        "set_point_PD2_AOM_A2",
        "set_point_PD3_AOM_A3",
        "set_point_PD4_AOM_A4",
        "set_point_PD5_AOM_A5",
        "set_point_PD6_AOM_A6",
    )

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        # Repository examination supplies None for arguments and may run with
        # an empty dataset namespace; bind the suffixed device superset now and
        # defer mode validation and dataset loading to prepare().
        self.base = BaseExperimentMasterSatellite(experiment=self)
        self.base.build()

        self.setattr_argument(
            "which_node",
            EnumerationValue(self.VALID_NODES),
            "Node selection",
            tooltip="Node1 = alice, Node2 = bob. The other node's beams, "
                    "coils and microwave/RF DDS are forced off for the run.",
        )

        ### overwrite the experiment variables of the same names
        self.setattr_argument("t_SPCM_exposure", NumberValue(10 * ms, unit='ms'))
        self.setattr_argument("atom_counts_per_s_threshold", NumberValue(14000))
        self.setattr_argument("target_cost", NumberValue(-3000.0))
        self.setattr_argument("n_measurements", NumberValue(10, type='int', scale=1, ndecimals=0, step=1))
        self.setattr_argument("set_best_parameters_at_finish", BooleanValue(True))
        self.both_mode = "coils and beam powers"
        self.beam_mode = "beam powers only"
        self.coil_mode = "coils only"
        self.setattr_argument("what_to_tune", EnumerationValue([self.both_mode, self.coil_mode, self.beam_mode]))

        group1 = "optimizer settings"
        self.setattr_argument("max_runs",NumberValue(100, type='int', scale=1, ndecimals=0, step=1),group1)
        self.setattr_argument(
            "compare_prev_and_optimized_values",
            BooleanValue(True),
            group1,
            tooltip="After the search, re-measure the cost at the previous "
                    "values and at the optimized values, and keep the "
                    "optimized ones only if they actually measured better. "
                    "Costs two extra measurement rounds.",
        )

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



        ### base.set_datasets_from_gui_args() has no master-satellite
        ### equivalent; the shadowing arguments are re-asserted in prepare().
        print("build - done")

    def prepare(self):
        ### t_SPCM_exposure and n_measurements are run-local GUI values that
        ### share their names with persistent variables (node-suffixed and
        ### global respectively). Capture them before configure_execution loads
        ### those datasets over the same attributes, and re-assert them after,
        ### so the submitted values win for this run.
        submitted_t_SPCM_exposure = self.t_SPCM_exposure
        submitted_n_measurements = self.n_measurements

        selected_node = str(self.which_node)
        if selected_node not in self.VALID_NODES:
            raise ValueError(
                f"Unsupported which_node {self.which_node!r}; expected "
                "'Node1' or 'Node2'."
            )
        self._selected_node = selected_node

        self.base.configure_execution("single_node", selected_node)
        self.base.prepare()

        ### The which_node ARGUMENT is consumed here and the attribute is
        ### replaced by the legacy presentation, exactly as
        ### MicrowaveScanOptimizer_master_satellite does it. The reused
        ### single-node code reads self.which_node expecting 'alice'/'bob':
        ### aom_feedback.py builds utilities/config/<which_node>/
        ### feedback_channels.json from it, and there is no config/Node1
        ### directory, so leaving 'Node1' here makes the stabilizer raise
        ### FileNotFoundError. self._selected_node keeps the Node1/Node2 form
        ### for anything that needs it.
        self.which_node = self.base.NODE_LEGACY_NAMES[selected_node]

        ### Must follow base.prepare(): it asserts the base is prepared, and it
        ### installs the node-suffixed feedback dataset map the stabilizer uses.
        ### Must also follow the which_node publication above, since the
        ### stabilizer reads it while being constructed.
        self.base.prepare_laser_stabilizer(
            stabilizer_factory=self._stabilizer_factory
        )

        self.t_SPCM_exposure = submitted_t_SPCM_exposure
        self.n_measurements = submitted_n_measurements

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

        ### Resolved to the selected node's authoritative names, because these
        ### are persisted at finish. resolve_experiment_variable_target is the
        ### resolver for experiment VARIABLES; the dataset redirect covers the
        ### coil volts but not the set points, so neither list may rely on it.
        self.volt_datasets = [
            self.base.resolve_experiment_variable_target(name)
            for name in self.COIL_VOLT_VARIABLES
        ]
        self.setpoint_datasets = [
            self.base.resolve_experiment_variable_target(name)
            for name in self.SET_POINT_VARIABLES
        ]
        ### Read back under the legacy projection, exactly as standalone does:
        ### in single_node mode the bare names are present on the experiment.
        self.default_setpoints = np.array(
            [getattr(self, name) for name in self.SET_POINT_VARIABLES]
        )

        self.atom_loading_time_list = np.zeros(self.n_measurements)
        self.set_dataset(self.AllSPCMs_rate_dataset,
                         [0.0],
                         broadcast=True)

        ### record_FORT_MM_power / record_FORT_APD_power (k10cr1_functions.py)
        ### APPEND to these without creating them, so the caller must seed
        ### them or the first append raises
        ###   KeyError: Cannot mutate nonexistent dataset 'FORT_MM_monitor_Node1'
        ### Seeded unconditionally and not "only if absent": append_to_dataset
        ### resolves through DatasetManager._get_mutation_target(), which sees
        ### only THIS worker's local and broadcaster dicts, while get() falls
        ### through to the master's database -- so a leftover dataset from an
        ### earlier run satisfies get() and still fails the append. Base says
        ### the same thing at initialize_single_node_result_state().
        ###
        ### Only these two, rather than calling
        ### base.initialize_single_node_result_state(): that method also sets
        ### atom_loading_time_list to a PYTHON LIST, and get_cost takes
        ### TArray(TFloat, 1), so it would undo the numpy array assigned just
        ### above. FORT_Polarization_Optimizer_master_satellite seeds exactly
        ### these two for the same reason.
        self.set_dataset("FORT_MM_monitor", [], broadcast=True)
        self.set_dataset("FORT_APD_monitor", [], broadcast=True)

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

        ### The settings that were in place BEFORE this optimization, expressed
        ### in M-LOOP's own parameterisation so the very same
        ### optimization_routine can re-measure them. Coil volts are absolute,
        ### so the originals go in as they are; set points are tuned as
        ### multipliers of default_setpoints, so "unchanged" is exactly 1.0.
        ### Built here, in bounds order, because optimization_routine
        ### overwrites self.coil_values with whatever it last tried -- by the
        ### end of the search the originals are no longer on the experiment.
        previous_params = []
        if self.tune_coils:
            previous_params += [self.AZ_bottom_volts_MOT,
                                self.AZ_top_volts_MOT,
                                self.AX_volts_MOT,
                                self.AY_volts_MOT]
        if self.tune_beams:
            n_tuned_beams = 4 if self.disable_z_beam_tuning else 6
            previous_params += [1.0] * n_tuned_beams
        self.previous_params = np.array(previous_params)

        print("max bounds")
        print(max_bounds)
        print("min bounds")
        print(min_bounds)

        self.mloop_controller = mlc.create_controller(interface,
                                           max_num_runs=self.max_runs,
                                           target_cost=self.target_cost, # -10000 corresponds to average atom_loading_time = 100ms calculated from -1000/atom_loading_time.
                                           num_params=n_params,
                                           min_boundary=min_bounds,
                                           max_boundary=max_bounds)
        print("prepare - done")

    def run(self):
        self.initialize_hardware()
        ### Leave the idle node dark: its beams, coils and microwave/RF DDS
        ### would otherwise keep whatever state the last experiment left, and
        ### its MOT light reaches this node's SPCMs.
        self.force_other_node_off()
        self.warm_up()

        self.mloop_controller.optimize()

        print('Best parameters found:')
        print(self.mloop_controller.best_params)
        best_params = self.mloop_controller.best_params

        ### Re-measure both settings before overwriting anything. M-LOOP's best
        ### cost belongs to whichever run was luckiest during the search, which
        ### is not a fair comparison against the settings already in place.
        if self.compare_prev_and_optimized_values:
            best_params = self.compare_previous_and_optimized(best_params)

        self.set_experiment_variables_to_best_params(best_params)

    @kernel
    def initialize_hardware(self):
        self.base.initialize_hardware()

    @kernel
    def force_other_node_off(self):
        self.base.force_other_node_off()

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

        ### warm up to make sure we get to the setpoints.
        ###
        ### break_realtime AFTER every run, which is the placement that
        ### matters. laser_stabilizer.run() advances the timeline only by its
        ### coded delays (0.1 ms either side of each measurement) while
        ### spending real time on sampler reads and host RPCs, and
        ### subroutines/aom_feedback.py manages no slack of its own. Measured
        ### on hardware, a single run RETURNS about 5.8 ms behind the wall
        ### clock, so the next RTIO event underflows whether that event is the
        ### following iteration or the FORT switch after the loop. Re-arming
        ### before the run only fixed the former: the first attempt underflowed
        ### inside the loop at -4.9 ms, the second on dds_FORT.sw.off() at
        ### -5.8 ms. Re-arming after covers both.
        ###
        ### Safe here precisely because warm-up measures nothing: it just
        ### drives the AOMs to their set points, so moving the cursor forward
        ### changes no interval that matters. Other master-satellite
        ### experiments get the same effect by calling the stabilizer once per
        ### kernel entry, right after core.reset().
        for i in range(10):
            self.laser_stabilizer.run()
            self.core.break_realtime()

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

            ### Same slack erosion as in warm_up, re-armed the same way: after
            ### the run, because the run returns with the cursor behind. The
            ### zotino write below is exactly the kind of RTIO event that
            ### underflowed otherwise. Every feedback call here happens BEFORE
            ### t_before_atom is read, so moving the cursor forward cannot
            ### affect a measured loading time.
            for i in range(3):
                self.laser_stabilizer.run()
                self.core.break_realtime()
        else:
            self.laser_stabilizer.run()
            self.core.break_realtime()

        if self.tune_coils:
            self.zotino0.set_dac(self.coil_values, channels=self.coil_channels)

        delay(1 * ms)

        ##################### This is the core of the optimizer that runs the sequence and get a cost:

        ### reset the counts dataset each run so we don't overwhelm the dashboard when plotting
        self.set_dataset(self.AllSPCMs_rate_dataset, [0.0], broadcast=True)

        ### The final feedback before measuring, and the one place this
        ### experiment records the FORT MM and APD powers. The standalone file
        ### calls bare laser_stabilizer.run() here, which feeds back to the 6
        ### MOT AOMs and the FORT LOADING setpoint only, so FORT_MM_monitor
        ### and FORT_APD_monitor were never written -- the GVS experiment
        ### functions get them precisely because they go through this helper.
        ###
        ### It also runs stabilizer_FORT at setpoints 2 and 1 (science holding
        ### and science) before the MOT feedback, so each M-LOOP run gains two
        ### extra stabilizer passes, and it leaves the FORT ON at the loading
        ### setpoint because it switches the FORT on to record. The
        ### measurement loop below switches the FORT on itself a few ms later,
        ### so the only difference is the FORT being on slightly earlier in
        ### the first measurement.
        ###
        ### Deliberately only here, not at the other feedback sites: warm_up
        ### runs feedback ten times, and this helper does three stabilizer
        ### passes instead of one, so putting it there would triple the
        ### warm-up cost for ten redundant MM points.
        run_feedback_and_record_FORT_MM_power(self)
        ### Redundant -- the helper already re-arms on return, deliberately, so
        ### "no caller can rely on a deterministic cursor across a feedback
        ### run". Kept so every feedback site in this file reads the same way.
        self.core.break_realtime()

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
                with parallel:
                    self.ttl_SPCM0_counter.gate_rising(atom_check_time)
                    self.ttl_SPCM1_counter.gate_rising(atom_check_time)
                    self.ttl_SPCM0_OtherNode_counter.gate_rising(atom_check_time)
                    self.ttl_SPCM1_OtherNode_counter.gate_rising(atom_check_time)

                AllSPCMs_atom_check = int(self.ttl_SPCM0_counter.fetch_count() + self.ttl_SPCM1_counter.fetch_count() + \
                                      self.ttl_SPCM0_OtherNode_counter.fetch_count() + self.ttl_SPCM1_OtherNode_counter.fetch_count())


                AllSPCMs_counts_per_s = AllSPCMs_atom_check / atom_check_time
                delay(1 * ms)
                self.append_to_dataset(self.AllSPCMs_rate_dataset, AllSPCMs_counts_per_s)
                try_n += 1

                if AllSPCMs_counts_per_s > self.atom_counts_per_s_threshold:
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

            ### Re-armed exactly as in optimization_routine; see warm_up.
            for i in range(3):
                self.laser_stabilizer.run()
                self.core.break_realtime()
        else:
            self.laser_stabilizer.run()
            self.core.break_realtime()

        # if self.tune_coils:
        #     self.zotino0.set_dac(self.coil_values, channels=self.coil_channels)

        delay(1 * ms)

        ##################### This is the core of the optimizer that runs the sequence and get a cost:

        ### reset the counts dataset each run so we don't overwhelm the dashboard when plotting
        self.set_dataset(self.AllSPCMs_rate_dataset, [0.0], broadcast=True)

        for i in range(self.n_measurements):
            if self.tune_coils:
                self.zotino0.set_dac(self.coil_values, channels=self.coil_channels)
            delay(10*ms)
            # delay(1 * ms)
            self.laser_stabilizer.run()
            ### The DDS switching below is the next RTIO event, and this is all
            ### before t_before_atom is read, so re-arming cannot affect the
            ### measured loading time.
            self.core.break_realtime()
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
                with parallel:
                    self.ttl_SPCM0_counter.gate_rising(atom_check_time)
                    self.ttl_SPCM1_counter.gate_rising(atom_check_time)
                    self.ttl_SPCM0_OtherNode_counter.gate_rising(atom_check_time)
                    self.ttl_SPCM1_OtherNode_counter.gate_rising(atom_check_time)

                AllSPCMs_atom_check = int(self.ttl_SPCM0_counter.fetch_count() + self.ttl_SPCM1_counter.fetch_count() + \
                                      self.ttl_SPCM0_OtherNode_counter.fetch_count() + self.ttl_SPCM1_OtherNode_counter.fetch_count())

                AllSPCMs_counts_per_s = AllSPCMs_atom_check / atom_check_time
                delay(1 * ms)
                self.append_to_dataset(self.AllSPCMs_rate_dataset, AllSPCMs_counts_per_s)
                try_n += 1

                if AllSPCMs_counts_per_s > self.atom_counts_per_s_threshold:
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


    def describe_parameters(self, params):
        """(authoritative name, absolute value) pairs for a parameter vector.

        Walks the vector in the same order prepare() built the bounds in, and
        converts set-point multipliers to the absolute values that would be
        persisted, so printed numbers match what lands in the datasets.
        """
        described = []
        index = 0
        if self.tune_coils:
            for name in self.volt_datasets:
                described.append((name, float(params[index])))
                index += 1
        if self.tune_beams:
            n_tuned_beams = 4 if self.disable_z_beam_tuning else 6
            for beam_i in range(n_tuned_beams):
                described.append((
                    self.setpoint_datasets[beam_i],
                    float(self.default_setpoints[beam_i] * params[index]),
                ))
                index += 1
        return described

    def compare_previous_and_optimized(self, optimized_params):
        """Re-measure the cost at the previous and the optimized settings.

        Returns the parameter vector that should actually be persisted. The
        optimized one is returned only if it clears BOTH bars:

          * strictly better than the previous values, re-measured here. A tie
            keeps the previous ones, since there is no reason to overwrite
            working settings for no gain.
          * better than target_cost. Beating the starting point is not the
            same as being good enough, and M-LOOP stops at max_runs whether
            or not the target was reached -- RID 38552 ended that way with
            -2793 against a target of -8000. Without this bar, a search that
            merely nudged the cost would overwrite tuned settings.

        A failure of the second bar is reported distinctly from a genuine
        non-improvement, because the remedy differs: more max_runs, or a
        target_cost that matches what the apparatus can actually reach.

        Note the bar lives here, in the comparison. With
        compare_prev_and_optimized_values unticked there is no re-measurement
        and no target check, and M-LOOP's best is persisted as the standalone
        file has always done.

        Both are measured with the same optimization_routine the search used,
        so the two costs are directly comparable: same sequence, same
        n_measurements, same averaging. Note the routine keeps its own
        best-cost bookkeeping, so these two runs also append to the cost
        dataset -- that is deliberate, it makes the comparison visible on the
        cost plot alongside the search.

        Order is previous first, then optimized. Anything that drifts slowly
        over the pair therefore counts against the optimized values, which is
        the conservative direction for a decision about whether to overwrite.
        """
        previous_params = self.previous_params

        if len(optimized_params) != len(previous_params):
            raise ValueError(
                f"The optimized parameter vector has {len(optimized_params)} "
                f"entries but the previous-value vector has "
                f"{len(previous_params)}. These are built from the same "
                "tuning mode and must match."
            )

        ### M-LOOP reports nan parameters when its learner fails. Applying
        ### those would drive the coils with nan and then persist it, so stop
        ### before the comparison rather than after.
        if not np.all(np.isfinite(optimized_params)):
            print("")
            print("COMPARISON SKIPPED: M-LOOP returned non-finite parameters "
                  f"{optimized_params}.")
            print("Keeping the previous values. Nothing is overwritten.")
            return previous_params

        print("")
        print("=" * 72)
        print("COMPARING the previous values against the optimized values")
        print(f"  Each is re-measured with the same sequence and the same "
              f"n_measurements = {self.n_measurements}.")
        print("  Cost is the average of -1000/t_load, so LOWER is better.")
        print("=" * 72)

        print("  re-measuring the PREVIOUS values ...")
        previous_cost = self.optimization_routine(previous_params)
        print(f"  previous  cost = {previous_cost}")

        print("  re-measuring the OPTIMIZED values ...")
        optimized_cost = self.optimization_routine(optimized_params)
        print(f"  optimized cost = {optimized_cost}")

        print("-" * 72)
        better_than_previous = optimized_cost < previous_cost
        meets_target = optimized_cost < self.target_cost

        if not better_than_previous:
            print(f"NOT IMPROVED: previous {previous_cost} vs optimized "
                  f"{optimized_cost}.")
            print("KEEPING THE PREVIOUS VALUES. Nothing is overwritten.")
            print("=" * 72)
            print("")
            return previous_params

        if not meets_target:
            print(f"BETTER BUT SHORT OF TARGET: cost {previous_cost} -> "
                  f"{optimized_cost}, which does not reach target_cost "
                  f"{self.target_cost:.0f}.")
            print("KEEPING THE PREVIOUS VALUES. Nothing is overwritten.")
            print("To accept a result like this, either give the search more "
                  "room with max_runs,")
            print("or relax target_cost to what the apparatus can actually "
                  "reach.")
            print("=" * 72)
            print("")
            return previous_params

        print(f"IMPROVED AND TARGET MET: cost {previous_cost} -> "
              f"{optimized_cost} (better by "
              f"{previous_cost - optimized_cost}, target "
              f"{self.target_cost:.0f}).")
        print("UPDATING to the optimized values:")
        for (name, old_value), (_, new_value) in zip(
                self.describe_parameters(previous_params),
                self.describe_parameters(optimized_params)):
            print(f"    {name:32} {old_value:12.6f} -> {new_value:12.6f}")
        print("=" * 72)
        print("")
        return optimized_params

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
        ### The M-LOOP plots are diagnostics, so a failure to draw them must
        ### not fail the run: by the time analyze() is reached the measurements
        ### are done and the decision about what to persist has been made.
        ###
        ### This is not hypothetical. M-LOOP 3.3.4 cannot plot a single-run
        ### optimization: when the target cost is reached on run 0 it stores
        ### in_costs as a 0-dimensional array, and plot_cost_vs_run then raises
        ### "IndexError: too many indices for array" (visualizations.py:469).
        ### RID 38545 completed its comparison correctly and was still marked
        ### failed for exactly that reason.
        try:
            mlv.show_all_default_visualizations(self.mloop_controller)
        except Exception as error:
            print(f"M-LOOP visualizations skipped: "
                  f"{type(error).__name__}: {error}")
            print("The optimization itself completed; only the plots failed.")
