"""Atom loading optimization for the master-satellite architecture.

Ported from standalone/AtomLoadingOptimizer_load_until_atom.py, which stays
untouched and keeps serving the standalone systems. MOT and FORT beams are
turned on and M-LOOP minimizes the time taken to load an atom, tuning the MOT
coil volts and/or the six MOT beam set points.

which_node selects Node1, Node2 or both_nodes. Node1/Node2 tune one node's
coils and beams against the SPCM rate with the other node forced dark.
both_nodes (2026-11-03) runs Base in two_nodes mode and tunes BOTH nodes
against one joint cost, loading the way the two-node experiments do; see BOTH
NODES at the end. (This used to say two_nodes was impossible because
master-satellite laser feedback was single-node only; Base has built one
stabilizer per node since e768d50.)

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
  the standalone file. Only one optimizer runs at a time, so there is nothing
  to collide with, and the applets read these names as they are.

BOTH NODES (which_node = both_nodes)

* Base is configured two_nodes, so no bare per-node device or variable names
  exist and the single-node kernels cannot compile. The mode has its own
  kernels, warm_up_both_nodes and optimization_routine_both_nodes, written out
  with explicit _Node1/_Node2 names, and the host picks which kernel to call
  (run, measure_cost). The single-node kernels are unchanged. Raw device names
  (zotino0, sampler1, ...) still exist bare in two_nodes mode and mean NODE1's
  crate, so they are never used here.

* Coils: the two-node coil set, <coil>_for_two_nodes_NodeX, is what is tuned
  and saved, whatever override_coils_with_two_nodes_setting says; that flag
  only decides whether the two-node experiments load with it (prepare warns
  when it is off). The single-node coils are never touched. When the coils
  are NOT tuned, the run loads with the coils the two-node experiments use
  now (Base's MOT_coil_volts_NodeX, chosen by the flag), so beam powers are
  tuned under the field they will be used with. Set points: each node's
  set_point_PD*_AOM_A*_NodeX
  (there is no separate two-node set). Parameter order: Node1's 4 coils,
  Node2's 4 coils, then the set-point multipliers Node1 A1..A4[,A5,A6] and
  Node2 A1..A4[,A5,A6]. The coil bounds must be differential (+-dV around each
  node's own values): the absolute bound arguments are one set and fit
  neither node.

* Each measurement loads the way load_until_atom_in_both_nodes_recycle
  (subroutines/experiment_functions_two_nodes.py) does: MOT coils written at
  the start of the load, both cooling DPs set to MOT, repump/cooling/A1-A6
  switched in groups of at most 4 events 5 us apart (a 9th event on one
  timestamp is silently dropped, RID 38676), both FORTs at the
  feedback-corrected loading amplitude, Node1's UV trigger, then 20 ms checks
  of the four SPCMs summed against atom_counts_per_s_threshold_both_nodes;
  afterwards UV off, MOT light off, coils to PGC. Not copied from the
  experiment: its rounds with feedback in between, and its reuse of surviving
  atoms. Here every measurement, the first included, starts from empty traps
  (both FORTs off) and a dissipated MOT (300 ms), and a load that has not
  passed after 100 checks scores as a failure (10e9 s), so the cost is the
  time to load both traps from empty.

* The joint check cannot tell which node holds an atom, so the cost follows
  the SLOWER node, and the faster node's parameters are only loosely pinned
  within +-dV. Every measurement also records its loading time and the summed
  rate at the first and at the last check (both_nodes_loading_time,
  both_nodes_first_check_counts_per_s, both_nodes_last_check_counts_per_s),
  so a cost that improved only because the background got brighter can be
  spotted afterwards.

* Feedback is experiment_functions_two_nodes.run_feedback_and_record_FORT_MM_
  power (both nodes; it has the single-node helper's name, so it is imported
  under an alias), recording into FORT_MM/APD_monitor_Node1/_Node2, which
  prepare seeds. Every laser_stabilizer_NodeX.run() is followed by
  core.break_realtime(), as in single-node mode.

* Saving: the comparison is joint, so both nodes' values are kept or reverted
  together. save_best_params_both_nodes writes them on the host, as the same
  name/value pairs describe_parameters prints.
"""

from artiq.experiment import *
import numpy as np
import logging

### Imports for M-LOOP. Only mloop.interfaces is imported here, because
### MLOOPInterface below subclasses it. mloop.controllers (it loads
### TensorFlow, ~8.6 s) is imported in prepare() and mloop.visualizations in
### analyze(), where they are used: ARTIQ 7 kills any build -- module import
### plus build() -- that takes over 15 s (artiq/master/worker.py,
### Worker.build timeout=15.0), and with all three here the build took ~10.6 s
### of M-LOOP imports alone. Submitted while another run kept the master busy,
### 40267 and 40268 (11-02) died in build. prepare() and analyze() have no
### such timeout.
import mloop.interfaces as mli

### Imported by name rather than with the `import *` MicrowaveScanOptimizer
### uses: this file needs exactly one function from that 9,500-line module, and
### naming it keeps the dependency visible.
from subroutines.experiment_functions import (
    run_feedback_and_record_FORT_MM_power,
)
### The two-node helper of the same name, for which_node = both_nodes: feedback
### on both nodes, FORT powers into FORT_MM/APD_monitor_Node1/_Node2.
from subroutines.experiment_functions_two_nodes import (
    run_feedback_and_record_FORT_MM_power
    as run_feedback_and_record_FORT_MM_power_both_nodes,
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

    Minimize atom loading time with M-LOOP, on one node or on both together
    (which_node = both_nodes), tuning MOT coil volts and/or MOT beam set
    points.
    """

    VALID_NODES = ("Node1", "Node2", "both_nodes")
    NODES = ("Node1", "Node2")
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
            tooltip="Node1 = alice, Node2 = bob: that node alone, with the "
                    "other node's beams, coils and microwave/RF DDS forced "
                    "off. both_nodes: both nodes load together, as in the "
                    "two-node experiments, and both nodes' parameters are "
                    "tuned (coils: the <coil>_for_two_nodes_NodeX set).",
        )

        ### overwrite the experiment variables of the same names
        self.setattr_argument("t_SPCM_exposure", NumberValue(10 * ms, unit='ms'))
        ### The loading check's threshold, counts/s summed over the four SPCMs,
        ### one per choice of which_node because the levels differ: with both
        ### MOTs on, the four SPCMs read ~55k counts/s with no atom and ~80k
        ### with one (upgrade_notes/two_atom_threshold_and_node1_fluorescence.md).
        self.setattr_argument("atom_counts_per_s_threshold_node1", NumberValue(55000))
        self.setattr_argument("atom_counts_per_s_threshold_node2", NumberValue(45000))
        self.setattr_argument("atom_counts_per_s_threshold_both_nodes", NumberValue(100000))
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
                    "target_cost is reported, not required. Costs two extra "
                    "measurement rounds.",
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
                "'Node1', 'Node2' or 'both_nodes'."
            )
        self._selected_node = selected_node
        self.both_nodes = selected_node == "both_nodes"

        if self.both_nodes:
            ### two_nodes mode takes no which_node. base.prepare() then sets
            ### self.which_node to 'alice' (Base._publish_broadcast_aliases);
            ### nothing here reads it, because each node's stabilizer is built
            ### with its own config_node. Branch on self.both_nodes instead.
            self.base.configure_execution("two_nodes", None)
            self.base.prepare()
        else:
            self.base.configure_execution("single_node", selected_node)
            self.base.prepare()

            ### The which_node ARGUMENT is consumed here and the attribute is
            ### replaced by the legacy presentation, exactly as
            ### MicrowaveScanOptimizer_master_satellite does it. The reused
            ### single-node code reads self.which_node expecting 'alice'/'bob':
            ### aom_feedback.py builds utilities/config/<which_node>/
            ### feedback_channels.json from it, and there is no config/Node1
            ### directory, so leaving 'Node1' here makes the stabilizer raise
            ### FileNotFoundError. self._selected_node keeps the Node1/Node2
            ### form for anything that needs it.
            self.which_node = self.base.NODE_LEGACY_NAMES[selected_node]

        ### Must follow base.prepare(): it asserts the base is prepared, and it
        ### installs the node-suffixed feedback dataset map the stabilizer uses.
        ### Must also follow the which_node publication above, since the
        ### stabilizer reads it while being constructed. In two_nodes mode it
        ### builds one stabilizer per node instead: laser_stabilizer_Node1/
        ### _Node2, stabilizer_AOM_A1_Node1, ..., stabilizer_FORT_Node2.
        self.base.prepare_laser_stabilizer(
            stabilizer_factory=self._stabilizer_factory
        )

        self.t_SPCM_exposure = submitted_t_SPCM_exposure
        self.n_measurements = submitted_n_measurements

        ### The threshold the kernels compare against, for this which_node.
        self.atom_counts_per_s_threshold = {
            "Node1": self.atom_counts_per_s_threshold_node1,
            "Node2": self.atom_counts_per_s_threshold_node2,
            "both_nodes": self.atom_counts_per_s_threshold_both_nodes,
        }[selected_node]
        self.atom_counts_threshold = self.atom_counts_per_s_threshold*self.t_SPCM_exposure

        if self.both_nodes:
            self._prepare_both_nodes()
            return

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
        ### the same thing at initialize_result_state().
        ###
        ### Only these two, rather than calling
        ### base.initialize_result_state(): that method also sets
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
        import mloop.controllers as mlc  # here, not at the top: see the imports
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

    def _check_saved_values_unchanged(self):
        """Refuse to start when a value this run would save changed after prepare.

        ARTIQ prepares a queued run while the run ahead of it is still
        running. An optimizer queued behind another one therefore read its
        starting values -- the centre of its coil bounds, the set points the
        beam search scales, and what "keep the previous values" writes back --
        before that run saved its result, and saving at the end would
        silently undo it. A both_nodes run saves both nodes' set points, so a
        Node1 or Node2 run followed by a both_nodes run is exactly that.
        Checked before anything touches the hardware: resubmit, and the new
        run starts from the new values.
        """
        changed = []
        for name, value in self.describe_parameters(self.previous_params):
            stored = float(self.get_dataset(name, archive=False))
            if stored != value:
                changed.append(f"{name} {value!r} at prepare, {stored!r} now")
        if changed:
            raise RuntimeError(
                "Values this run would save changed after it was prepared, "
                "most likely saved by a run queued ahead of it: "
                + "; ".join(changed) + ". Nothing was run or written; "
                "resubmit to start from the new values.")

    def _prepare_both_nodes(self):
        """The rest of prepare() for which_node = both_nodes.

        Kept apart from the single-node lines, which read the bare per-node
        names a two_nodes Base does not publish. Same steps in the same
        order: saved names and starting values, dataset seeding, bounds,
        M-LOOP. See BOTH NODES in the module docstring.
        """
        nodes = self.NODES
        if not self.use_differential_boundaries:
            raise ValueError(
                "which_node = both_nodes needs use_differential_boundaries: "
                "the absolute coil bounds (V_AZ_bottom_min, ...) are one set "
                "for both nodes, and the two nodes' MOT coils are volts apart "
                "(AZ_bottom ~6.5 V on Node1, ~4.7 V on Node2 on 2026-11-03)."
            )

        self.tune_beams = self.what_to_tune == self.beam_mode or self.what_to_tune == self.both_mode
        self.tune_coils = self.what_to_tune == self.coil_mode or self.what_to_tune == self.both_mode
        self.n_tuned_beams_per_node = 4 if self.disable_z_beam_tuning else 6

        ### The two-node coil set: tuned and saved whatever
        ### override_coils_with_two_nodes_setting says (the flag only decides
        ### whether the two-node experiments load with it), so the
        ### single-node coils are never touched.
        self.volt_datasets = [
            self.base.resolve_experiment_variable_target(
                f"{name}_for_two_nodes", node)
            for node in nodes
            for name in self.COIL_VOLT_VARIABLES
        ]
        start_volts = [float(getattr(self, name)) for name in self.volt_datasets]
        if self.tune_coils:
            loading_volts = start_volts
            if not getattr(self, self.base.TWO_NODE_MOT_COIL_FLAG):
                print("WARNING both_nodes: tuning the two-node coil set "
                      "(<coil>_for_two_nodes_NodeX), but "
                      f"{self.base.TWO_NODE_MOT_COIL_FLAG} is False, so the "
                      "two-node experiments keep loading with each node's "
                      "single-node coils until it is switched on.")
        else:
            ### Not tuned: load with the coils the two-node experiments use
            ### now -- the single-node set, or the two-node set when the flag
            ### is on (Base) -- so the beam powers are tuned under that field.
            loading_volts = (list(self.MOT_coil_volts_Node1)
                             + list(self.MOT_coil_volts_Node2))
        ### Overwritten by optimization_routine_both_nodes with each tried value.
        self.coil_values_Node1 = np.array(loading_volts[:4])
        self.coil_values_Node2 = np.array(loading_volts[4:])
        self.setpoint_datasets = [
            self.base.resolve_experiment_variable_target(name, node)
            for node in nodes
            for name in self.SET_POINT_VARIABLES
        ]
        self.default_setpoints = np.array(
            [float(getattr(self, name)) for name in self.setpoint_datasets]
        )

        self.atom_loading_time_list = np.zeros(self.n_measurements)
        self.set_dataset(self.AllSPCMs_rate_dataset, [0.0], broadcast=True)
        ### record_FORT_powers_both_nodes APPENDS to these without creating
        ### them, and in two_nodes mode the redirect suffixes nothing, so they
        ### are seeded under their per-node names. (Base.initialize_result_
        ### state seeds them for the GVS, but it would also turn
        ### atom_loading_time_list into a Python list.)
        for node in nodes:
            self.set_dataset(f"FORT_MM_monitor_{node}", [], broadcast=True)
            self.set_dataset(f"FORT_APD_monitor_{node}", [], broadcast=True)
        ### One entry per measurement over the whole run.
        for name in ("both_nodes_loading_time",
                     "both_nodes_first_check_counts_per_s",
                     "both_nodes_last_check_counts_per_s"):
            self.set_dataset(name, [], broadcast=True)

        self.cost_dataset = "cost"
        self.current_best_cost = 0
        self.set_dataset(self.cost_dataset,
                         [self.current_best_cost],
                         broadcast=True)

        ### Bounds and the pre-optimization values, in parameter order:
        ### Node1's 4 coils, Node2's 4 coils, then the set-point multipliers
        ### Node1 A1..A4[,A5,A6] and Node2 A1..A4[,A5,A6]. Coils are absolute;
        ### set points are multipliers of default_setpoints, so "unchanged" is
        ### 1.0. Built here because the search overwrites coil_values_NodeX.
        min_bounds = []
        max_bounds = []
        previous_params = []
        if self.tune_coils:
            print("MLOOP will tune both nodes' coil volts (the _for_two_nodes set)")
            dV = [self.dV_AZ_bottom, self.dV_AZ_top, self.dV_AX, self.dV_AY] * 2
            min_bounds += [volts - d for volts, d in zip(start_volts, dV)]
            max_bounds += [volts + d for volts, d in zip(start_volts, dV)]
            previous_params += start_volts
        if self.tune_beams:
            print("MLOOP will tune both nodes' beam powers")
            n_beams = 2 * self.n_tuned_beams_per_node
            min_bounds += [1 - self.max_set_point_percent_deviation_minus] * n_beams
            max_bounds += [1 + self.max_set_point_percent_deviation_plus] * n_beams
            previous_params += [1.0] * n_beams
        n_params = len(max_bounds)
        self.best_params = np.zeros(n_params)
        self.previous_params = np.array(previous_params)

        print(f"both_nodes: {n_params} parameters; M-LOOP's learner starts "
              f"after max(10, 2 x {n_params}) = {max(10, 2 * n_params)} "
              f"training runs, of max_runs = {self.max_runs}.")
        print(f"both_nodes: a load passes at "
              f"{self.atom_counts_per_s_threshold:.0f} counts/s summed over "
              f"the four SPCMs (the two-node experiments use "
              f"two_atom_threshold_for_loading = "
              f"{self.two_atom_threshold_for_loading:.0f}).")
        if self.target_cost > -5000:
            print(f"WARNING both_nodes: target_cost is {self.target_cost:.0f}. "
                  "M-LOOP ends the search at the first evaluation at or below "
                  "target_cost, and the joint cost is noisy (about -1900 +- "
                  "525 at n_measurements 10, from RID 40342), so -3000 is "
                  "reached by chance in ~3% of evaluations. Use target_cost "
                  "<= -5000 and n_measurements >= 20.")

        import mloop.controllers as mlc  # here, not at the top: see the imports
        interface = MLOOPInterface()
        interface.get_next_cost_dict = self.get_next_cost_dict_for_mloop

        print("max bounds")
        print(max_bounds)
        print("min bounds")
        print(min_bounds)

        self.mloop_controller = mlc.create_controller(interface,
                                           max_num_runs=self.max_runs,
                                           target_cost=self.target_cost,
                                           num_params=n_params,
                                           min_boundary=min_bounds,
                                           max_boundary=max_bounds)
        print("prepare - done")

    def run(self):
        if self.set_best_parameters_at_finish:
            self._check_saved_values_unchanged()
        self.initialize_hardware()
        if self.both_nodes:
            ### initialize_hardware initialized both crates (every DDS off,
            ### both Zotinos at 0 V); there is no idle node to force off.
            self.warm_up_both_nodes()
        else:
            ### Leave the idle node dark: its beams, coils and microwave/RF
            ### DDS would otherwise keep whatever state the last experiment
            ### left, and its MOT light reaches this node's SPCMs.
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

        if self.both_nodes:
            self.save_best_params_both_nodes(best_params)
        else:
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

    @kernel
    def warm_up_both_nodes(self):
        """warm_up for which_node = both_nodes: both nodes' MOT light and
        FORTs on with the two-node MOT coils, ten feedback passes per node,
        then the light off again. Switch events in groups of at most 4, 5 us
        apart, as in load_until_atom_in_both_nodes_recycle."""

        self.core.reset()

        delay(1*ms)

        ### Set the coils to the two-node MOT loading setting
        self.zotino0_Node1.set_dac(
            [self.coil_values_Node1[0], self.coil_values_Node1[1],
             self.coil_values_Node1[2], self.coil_values_Node1[3]],
            channels=self.coil_channels_Node1)
        self.zotino0_Node2.set_dac(
            [self.coil_values_Node2[0], self.coil_values_Node2[1],
             self.coil_values_Node2[2], self.coil_values_Node2[3]],
            channels=self.coil_channels_Node2)
        delay(1 * ms)

        ### set the cooling DP AOMs to the MOT settings
        self.dds_cooling_DP_Node1.set(frequency=self.f_cooling_DP_MOT_Node1,
                                      amplitude=self.ampl_cooling_DP_MOT_Node1)
        self.dds_cooling_DP_Node2.set(frequency=self.f_cooling_DP_MOT_Node2,
                                      amplitude=self.ampl_cooling_DP_MOT_Node2)

        self.ttl_repump_switch_Node1.off()  ### turn on MOT RP
        self.ttl_repump_switch_Node2.off()
        delay(5 * us)
        self.dds_cooling_DP_Node1.sw.on()  ### turn on cooling
        self.dds_cooling_DP_Node2.sw.on()
        self.dds_AOM_A1_Node1.sw.on()
        self.dds_AOM_A1_Node2.sw.on()
        delay(5 * us)
        self.dds_AOM_A2_Node1.sw.on()
        self.dds_AOM_A2_Node2.sw.on()
        self.dds_AOM_A3_Node1.sw.on()
        self.dds_AOM_A3_Node2.sw.on()
        delay(5 * us)
        self.dds_AOM_A4_Node1.sw.on()
        self.dds_AOM_A4_Node2.sw.on()
        self.dds_AOM_A5_Node1.sw.on()
        self.dds_AOM_A5_Node2.sw.on()
        delay(5 * us)
        self.dds_AOM_A6_Node1.sw.on()
        self.dds_AOM_A6_Node2.sw.on()
        delay(0.1 * ms)

        self.dds_FORT_Node1.set(frequency=self.f_FORT_Node1,
                                amplitude=self.stabilizer_FORT_Node1.amplitude)
        self.dds_FORT_Node2.set(frequency=self.f_FORT_Node2,
                                amplitude=self.stabilizer_FORT_Node2.amplitude)
        self.dds_FORT_Node1.sw.on()
        self.dds_FORT_Node2.sw.on()

        ### delay for AOMs to thermalize
        delay(500 * ms)

        ### warm up to make sure we get to the setpoints; break_realtime
        ### AFTER every run, for the reason given in warm_up.
        for i in range(10):
            self.laser_stabilizer_Node1.run()
            self.core.break_realtime()
            self.laser_stabilizer_Node2.run()
            self.core.break_realtime()

        ### Turning off AOMs to be ready to start atom loading from scratch
        self.dds_FORT_Node1.sw.off()  ### turn off FORT
        self.dds_FORT_Node2.sw.off()
        self.dds_cooling_DP_Node1.sw.off()  ### turn off cooling
        self.dds_cooling_DP_Node2.sw.off()
        delay(5 * us)
        self.ttl_repump_switch_Node1.on()  ### turn off MOT RP
        self.ttl_repump_switch_Node2.on()


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
        ### Re-armed because of the set_dataset just above: from a kernel it
        ### is a BLOCKING RPC here (_DatasetRedirectMixin.set_dataset has no
        ### async flag, unlike stock ARTIQ's), so the cursor stands still
        ### while the host and the master answer, and the delay(1 ms) above
        ### was its only margin. 11 runs between 10-02 and 10-31 died between
        ### here and the end of the helper's first stabilizer pass; the one
        ### traceback kept is an RTIOUnderflow on its first event, the FORT
        ### dds.set (channel 36, slack -2432 mu). The 0.1 ms is the margin
        ### the GVS call sites give this helper. Still before t_before_atom,
        ### so no measured loading time changes.
        self.core.break_realtime()
        delay(0.1 * ms)

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
                param_idx = 4  # the set points follow the 4 coil values
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
    def optimization_routine_both_nodes(self, params: TArray(TFloat)) -> TInt32:
        """
        optimization_routine for which_node = both_nodes; see BOTH NODES in
        the module docstring.

        params: Node1's 4 coil volts and Node2's 4 (when tuning coils), then
            the set-point multipliers Node1 A1..A4[,A5,A6] and Node2
            A1..A4[,A5,A6] (when tuning beams).
        return:
            cost: the cost for the optimizer
        """

        self.core.reset()
        delay(1*ms)

        n_beams = self.n_tuned_beams_per_node
        if self.tune_coils:
            self.coil_values_Node1 = params[:4]
            self.coil_values_Node2 = params[4:8]
            if self.tune_beams:
                setpoint_multipliers = params[8:]
            else:
                setpoint_multipliers = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                                 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
        else:
            setpoint_multipliers = params

        if self.tune_beams:
            ### default_setpoints: Node1 A1..A6 at 0-5, Node2 A1..A6 at 6-11
            self.stabilizer_AOM_A1_Node1.set_points[0] = self.default_setpoints[0] * setpoint_multipliers[0]
            self.stabilizer_AOM_A2_Node1.set_points[0] = self.default_setpoints[1] * setpoint_multipliers[1]
            self.stabilizer_AOM_A3_Node1.set_points[0] = self.default_setpoints[2] * setpoint_multipliers[2]
            self.stabilizer_AOM_A4_Node1.set_points[0] = self.default_setpoints[3] * setpoint_multipliers[3]
            self.stabilizer_AOM_A1_Node2.set_points[0] = self.default_setpoints[6] * setpoint_multipliers[n_beams]
            self.stabilizer_AOM_A2_Node2.set_points[0] = self.default_setpoints[7] * setpoint_multipliers[n_beams + 1]
            self.stabilizer_AOM_A3_Node2.set_points[0] = self.default_setpoints[8] * setpoint_multipliers[n_beams + 2]
            self.stabilizer_AOM_A4_Node2.set_points[0] = self.default_setpoints[9] * setpoint_multipliers[n_beams + 3]
            if not self.disable_z_beam_tuning:
                self.stabilizer_AOM_A5_Node1.set_points[0] = self.default_setpoints[4] * setpoint_multipliers[4]
                self.stabilizer_AOM_A6_Node1.set_points[0] = self.default_setpoints[5] * setpoint_multipliers[5]
                self.stabilizer_AOM_A5_Node2.set_points[0] = self.default_setpoints[10] * setpoint_multipliers[n_beams + 4]
                self.stabilizer_AOM_A6_Node2.set_points[0] = self.default_setpoints[11] * setpoint_multipliers[n_beams + 5]

            ### Re-armed after every run, as in optimization_routine.
            for i in range(3):
                self.laser_stabilizer_Node1.run()
                self.core.break_realtime()
                self.laser_stabilizer_Node2.run()
                self.core.break_realtime()
        else:
            self.laser_stabilizer_Node1.run()
            self.core.break_realtime()
            self.laser_stabilizer_Node2.run()
            self.core.break_realtime()

        delay(1 * ms)

        ##################### This is the core of the optimizer that runs the sequence and get a cost:

        ### reset the counts dataset each run so we don't overwhelm the dashboard when plotting
        self.set_dataset(self.AllSPCMs_rate_dataset, [0.0], broadcast=True)
        ### Re-armed after this blocking RPC, as in optimization_routine.
        self.core.break_realtime()
        delay(0.1 * ms)

        ### The final feedback, both nodes, and the FORT MM and APD powers of
        ### both nodes. It leaves the MOT light and both FORTs on.
        run_feedback_and_record_FORT_MM_power_both_nodes(self)
        ### Redundant -- the helper re-arms on return -- kept so every feedback
        ### site in this file reads the same way.
        self.core.break_realtime()

        ### Every measurement, the first included, starts from empty traps and
        ### a dissipated MOT: the feedback above leaves the MOT light and the
        ### FORTs on, which would give measurement 0 a head start.
        self.ttl_repump_switch_Node1.on()  ### turn off MOT RP
        self.ttl_repump_switch_Node2.on()
        self.dds_cooling_DP_Node1.sw.off()  ### turn off cooling
        self.dds_cooling_DP_Node2.sw.off()
        delay(5 * us)
        self.dds_FORT_Node1.sw.off()  ### turn off FORT
        self.dds_FORT_Node2.sw.off()
        delay(1 * ms)
        ### set the coils to PGC settings, as the loader leaves them
        self.zotino0_Node1.set_dac(
            [self.AZ_bottom_volts_PGC_Node1, -self.AZ_bottom_volts_PGC_Node1,
             self.AX_volts_PGC_Node1, self.AY_volts_PGC_Node1],
            channels=self.coil_channels_Node1)
        self.zotino0_Node2.set_dac(
            [self.AZ_bottom_volts_PGC_Node2, -self.AZ_bottom_volts_PGC_Node2,
             self.AX_volts_PGC_Node2, self.AY_volts_PGC_Node2],
            channels=self.coil_channels_Node2)
        delay(300 * ms)  ### to dissipate the MOTs

        for i in range(self.n_measurements):
            ### Set the coils to the two-node MOT loading setting, at the start
            ### of every load as the loader does
            self.zotino0_Node1.set_dac(
                [self.coil_values_Node1[0], self.coil_values_Node1[1],
                 self.coil_values_Node1[2], self.coil_values_Node1[3]],
                channels=self.coil_channels_Node1)
            self.zotino0_Node2.set_dac(
                [self.coil_values_Node2[0], self.coil_values_Node2[1],
                 self.coil_values_Node2[2], self.coil_values_Node2[3]],
                channels=self.coil_channels_Node2)
            delay(1 * ms)

            ### set the cooling DP AOMs to the MOT settings
            self.dds_cooling_DP_Node1.set(frequency=self.f_cooling_DP_MOT_Node1,
                                          amplitude=self.ampl_cooling_DP_MOT_Node1)
            self.dds_cooling_DP_Node2.set(frequency=self.f_cooling_DP_MOT_Node2,
                                          amplitude=self.ampl_cooling_DP_MOT_Node2)

            self.ttl_repump_switch_Node1.off()  ### turn on MOT RP
            self.ttl_repump_switch_Node2.off()
            delay(5 * us)
            self.dds_cooling_DP_Node1.sw.on()  ### turn on cooling
            self.dds_cooling_DP_Node2.sw.on()
            self.dds_AOM_A1_Node1.sw.on()
            self.dds_AOM_A1_Node2.sw.on()
            delay(5 * us)
            self.dds_AOM_A2_Node1.sw.on()
            self.dds_AOM_A2_Node2.sw.on()
            self.dds_AOM_A3_Node1.sw.on()
            self.dds_AOM_A3_Node2.sw.on()
            delay(5 * us)
            self.dds_AOM_A4_Node1.sw.on()
            self.dds_AOM_A4_Node2.sw.on()
            self.dds_AOM_A5_Node1.sw.on()
            self.dds_AOM_A5_Node2.sw.on()
            delay(5 * us)
            self.dds_AOM_A6_Node1.sw.on()
            self.dds_AOM_A6_Node2.sw.on()
            delay(0.1 * ms)

            ### turn on the FORTs at their loading setpoints
            self.dds_FORT_Node1.set(frequency=self.f_FORT_Node1,
                                    amplitude=self.stabilizer_FORT_Node1.amplitude)
            self.dds_FORT_Node2.set(frequency=self.f_FORT_Node2,
                                    amplitude=self.stabilizer_FORT_Node2.amplitude)
            self.dds_FORT_Node1.sw.on()
            self.dds_FORT_Node2.sw.on()
            delay(5 * us)

            ### UV trigger is Node1 only; there is no Node2 UV pulse
            self.zotino0_Node1.set_dac([3.5], self.UV_trig_channel_Node1)
            delay(1 * ms)

            max_tries = 100  ### Maximum number of attempts before giving up
            atom_check_time = 20 * ms
            atom_loaded = False
            try_n = 0
            first_check_counts_per_s = 0.0
            last_check_counts_per_s = 0.0
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
                if try_n == 0:
                    first_check_counts_per_s = AllSPCMs_counts_per_s
                last_check_counts_per_s = AllSPCMs_counts_per_s
                try_n += 1

                if AllSPCMs_counts_per_s > self.atom_counts_per_s_threshold:
                    delay(100 * us)  ### Needs a delay of about 100us or maybe less
                    atom_loaded = True

            if atom_loaded:
                t_after_atom = now_mu()
                atom_loading_time = self.core.mu_to_seconds(t_after_atom - t_before_atom)
            else:
                atom_loading_time = 10e9 ### Just a large number to show no atom loading. This is compared to the typical 0.5 second atom loading.

            ### UV off, MOTs away
            self.zotino0_Node1.set_dac([0.0], self.UV_trig_channel_Node1)
            delay(1 * ms)
            self.atom_loading_time_list[i] = atom_loading_time
            self.append_to_dataset("both_nodes_loading_time", atom_loading_time)
            self.append_to_dataset("both_nodes_first_check_counts_per_s", first_check_counts_per_s)
            self.append_to_dataset("both_nodes_last_check_counts_per_s", last_check_counts_per_s)

            self.ttl_repump_switch_Node1.on()  ### turn off MOT RP
            self.ttl_repump_switch_Node2.on()
            self.dds_cooling_DP_Node1.sw.off()  ### turn off cooling
            self.dds_cooling_DP_Node2.sw.off()
            delay(5 * us)
            ### drop both atoms: the next load starts from empty traps
            self.dds_FORT_Node1.sw.off()  ### turn off FORT
            self.dds_FORT_Node2.sw.off()
            delay(1 * ms)
            ### set the coils to PGC settings, as the loader leaves them
            self.zotino0_Node1.set_dac(
                [self.AZ_bottom_volts_PGC_Node1, -self.AZ_bottom_volts_PGC_Node1,
                 self.AX_volts_PGC_Node1, self.AY_volts_PGC_Node1],
                channels=self.coil_channels_Node1)
            self.zotino0_Node2.set_dac(
                [self.AZ_bottom_volts_PGC_Node2, -self.AZ_bottom_volts_PGC_Node2,
                 self.AX_volts_PGC_Node2, self.AY_volts_PGC_Node2],
                channels=self.coil_channels_Node2)
            delay(300 * ms)  ### to dissipate the MOTs

        cost = self.get_cost(self.atom_loading_time_list)
        self.append_to_dataset(self.cost_dataset, cost)

        ################################################################################

        if cost < self.current_best_cost:
            self.current_best_cost = cost
            self.print_async("NEW BEST COST:", cost)
            param_idx = 0
            if self.tune_coils:
                self.print_async("BEST coil values:", params[:8])
                self.best_params[:8] = params[:8]
                param_idx = 8  # the set points follow the 8 coil values
            if self.tune_beams:
                for i in range(n_beams):
                    self.best_params[param_idx + i] = self.default_setpoints[i] * setpoint_multipliers[i]
                    self.best_params[param_idx + n_beams + i] = self.default_setpoints[6 + i] * setpoint_multipliers[n_beams + i]
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
        ### Re-armed exactly as after the same set_dataset in
        ### optimization_routine: with tune_coils the zotino write below is
        ### the next RTIO event, behind the same 1 ms margin.
        self.core.break_realtime()
        delay(0.1 * ms)

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
                param_idx = 4  # the set points follow the 4 coil values
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
            ### setpoint_datasets and default_setpoints hold A1..A6 of each
            ### node in turn: one node, or Node1 then Node2 for both_nodes.
            for node_start in range(0, len(self.setpoint_datasets), 6):
                for beam_i in range(node_start, node_start + n_tuned_beams):
                    described.append((
                        self.setpoint_datasets[beam_i],
                        float(self.default_setpoints[beam_i] * params[index]),
                    ))
                    index += 1
        return described

    def compare_previous_and_optimized(self, optimized_params):
        """Re-measure the cost at the previous and the optimized settings.

        Returns the parameter vector that should actually be persisted. The
        optimized one is returned when it is strictly better than the previous
        values, re-measured here. A tie keeps the previous ones, since there is
        no reason to overwrite working settings for no gain.

        target_cost is REPORTED, not required (since 2026-10-16, at the user's
        request). M-LOOP stops at max_runs whether or not the target was
        reached -- RID 38552 ended with -2793 against a target of -8000 --
        and an improvement over the settings in place is worth keeping even
        then. It used to be a second bar, which threw such improvements away.
        A result short of target is applied and flagged, with the remedy for
        getting further: more max_runs, or a target_cost the apparatus can
        actually reach.

        With compare_prev_and_optimized_values unticked there is no
        re-measurement at all, and M-LOOP's best is persisted as the
        standalone file has always done.

        Both are measured with the same cost routine the search used
        (measure_cost), so the two costs are directly comparable: same sequence, same
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
        previous_cost = self.measure_cost(previous_params)
        print(f"  previous  cost = {previous_cost}")

        print("  re-measuring the OPTIMIZED values ...")
        optimized_cost = self.measure_cost(optimized_params)
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

        if meets_target:
            print(f"IMPROVED AND TARGET MET: cost {previous_cost} -> "
                  f"{optimized_cost} (better by "
                  f"{previous_cost - optimized_cost}, target "
                  f"{self.target_cost:.0f}).")
        else:
            ### Better than what is in place, so it is applied even though the
            ### target was not reached; the target is reported, not enforced.
            print(f"IMPROVED BUT SHORT OF TARGET: cost {previous_cost} -> "
                  f"{optimized_cost} (better by "
                  f"{previous_cost - optimized_cost}), target_cost "
                  f"{self.target_cost:.0f} not reached.")
            print("Applying it anyway, since it beats the previous values. To "
                  "get further, give the search more max_runs or set a "
                  "target_cost the apparatus can reach.")
        print("UPDATING to the optimized values:")
        for (name, old_value), (_, new_value) in zip(
                self.describe_parameters(previous_params),
                self.describe_parameters(optimized_params)):
            print(f"    {name:32} {old_value:12.6f} -> {new_value:12.6f}")
        print("=" * 72)
        print("")
        return optimized_params

    def measure_cost(self, params):
        """One cost evaluation, with the kernel for this run's which_node."""
        if self.both_nodes:
            return self.optimization_routine_both_nodes(params)
        return self.optimization_routine(params)

    def get_next_cost_dict_for_mloop(self,params_dict):

        ### Get parameters from the provided dictionary
        params = params_dict['params']

        cost = self.measure_cost(params)
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

    def save_best_params_both_nodes(self, best_params):
        """set_experiment_variables_to_best_params for which_node = both_nodes.

        On the host: there is nothing to time, and describe_parameters already
        gives each value under the name it belongs to -- the two-node coil
        set and each node's set points -- so the table printed by the
        comparison and what is written are the same pairs.
        """
        if not self.set_best_parameters_at_finish:
            print("set_best_parameters_at_finish is off: nothing is saved.")
            return
        for name, value in self.describe_parameters(best_params):
            self.set_dataset(name, value, broadcast=True, persist=True)
            print(f"saved {name} = {value:.6f}")

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
            import mloop.visualizations as mlv  # see the imports at the top
            mlv.show_all_default_visualizations(self.mloop_controller)
        except Exception as error:
            print(f"M-LOOP visualizations skipped: "
                  f"{type(error).__name__}: {error}")
            print("The optimization itself completed; only the plots failed.")
