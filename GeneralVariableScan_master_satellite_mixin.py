"""Shared scan implementation (mixins) for the master-satellite GVS files.

This module holds no public EnvExperiment class and is not an ARTIQ Explorer
experiment; the runnable experiments live in the mode-named GVS files.
"""

import inspect
import logging
import re
import time

import numpy as np
from numpy import array

from artiq.experiment import *
from artiq.coredevice.exceptions import RTIOUnderflow

from utilities.BaseExperiment_master_satellite import (
    BaseExperimentMasterSatellite,
    merge_per_node_override_dictionaries,
    _DatasetRedirectMixin,
)


HISTORICAL_INDEPENDENT_TWO_NODE_EXPERIMENTS = frozenset({
    "Two_nodes_atom_loading_experiment",
    "Two_node_single_photon_experiment",
    "Two_node_single_photon_2_experiment",
    "Two_node_single_photon_2_optimization_experiment",
})


def build_experiment_function_registry(module, excluded_names=()):
    """Return locally defined experiment functions from *module*."""
    excluded_names = frozenset(excluded_names)
    return {
        name: function
        for name, function in vars(module).items()
        if inspect.isfunction(function)
        and function.__module__ == module.__name__
        and "experiment" in name
        and name not in excluded_names
    }


def build_single_node_function_registry():
    import subroutines.experiment_functions as experiment_functions

    return build_experiment_function_registry(
        experiment_functions,
        HISTORICAL_INDEPENDENT_TWO_NODE_EXPERIMENTS,
    )


def build_two_node_function_registry():
    import subroutines.experiment_functions_two_nodes as functions

    return build_experiment_function_registry(functions)


def snapshot_list_datasets(dataset_mgr):
    """Every list dataset this run holds, and its length: {key: (list, length)}.

    *dataset_mgr* is the run's ARTIQ DatasetManager. Archived datasets live in
    its `local`, broadcast ones in its broadcaster's raw_view (the same object
    when both). Taken before each attempt at a scan point, so that
    roll_back_list_datasets can take out what a failed attempt appended.
    """
    snapshot = {}
    for key, value in dataset_mgr.local.items():
        if isinstance(value, list):
            snapshot[key] = (value, len(value))
    for key, (_, value) in dataset_mgr._broadcaster.raw_view.items():
        if isinstance(value, list) and key not in snapshot:
            snapshot[key] = (value, len(value))
    return snapshot


def roll_back_list_datasets(dataset_mgr, snapshot):
    """Cut each list in *snapshot* that has grown since back to its length.

    Returns {key: entries removed}.

    Per-measurement series grow in place, through append_to_dataset, so a
    list that is still the SAME object and longer holds the failed attempt's
    entries at its end. A list replaced since (set_dataset, e.g. the
    _current_iteration histograms at every scan point) is left alone, and so
    is a dataset created since: the retry sets those again itself.

    The cut list is set again with the flags it has, so the master's dataset
    DB -- and with it the dashboard -- gets it too. A broadcast dataset must
    stay broadcast: a local write to it would delete it from the master DB.
    """
    broadcast_view = dataset_mgr._broadcaster.raw_view
    removed = {}
    for key, (value, length) in snapshot.items():
        if key in dataset_mgr.local:
            current = dataset_mgr.local[key]
        elif key in broadcast_view:
            current = broadcast_view[key][1]
        else:
            continue
        if current is not value or len(value) <= length:
            continue
        removed[key] = len(value) - length
        del value[length:]
        broadcast = key in broadcast_view
        dataset_mgr.set(
            key, value,
            broadcast=broadcast,
            persist=broadcast and broadcast_view[key][0],
            archive=key in dataset_mgr.local,
        )
    return removed


def _describe_rollback(removed):
    """'none', or e.g. '37 entries from AllSPCMs_RO1 and 41 other datasets'."""
    if not removed:
        return "none"
    key = max(sorted(removed), key=removed.get)
    others = len(removed) - 1
    return (f"{removed[key]} entries from {key}"
            + (f" and {others} other dataset(s)" if others else ""))


# ARTIQ 7 RTIOUnderflow texts: a CPU write (firmware ksupport/rtio.rs) reports
# the slack, a DMA playback (ksupport/lib.rs, dma_playback) does not.
_UNDERFLOW_MESSAGE = re.compile(
    r"RTIO underflow at (-?\d+) mu, channel (\d+)(?:, slack (-?\d+) mu)?")
_CORE_TRACEBACK_FRAME = re.compile(
    r'File "([^"]+)", line (\d+)(?:, column \d+)?, in ([^\s(]+)')


def parse_underflow_message(text):
    """(timestamp_mu, channel, slack_mu) from an RTIOUnderflow message, with
    slack_mu None for a DMA playback; None if the text does not parse."""
    match = _UNDERFLOW_MESSAGE.search(text)
    if match is None:
        return None
    timestamp, channel, slack = match.groups()
    return int(timestamp), int(channel), None if slack is None else int(slack)


def repository_frames(core_traceback):
    """'file:line function' of each frame from this repository in a core
    device traceback (str of the exception's artiq_core_exception), in the
    order ARTIQ prints them."""
    frames = []
    for path, line, function in _CORE_TRACEBACK_FRAME.findall(core_traceback):
        path = path.replace("\\", "/")
        if "qn_artiq_routines/" in path:
            frames.append(f"{path.rsplit('/', 1)[-1]}:{line} "
                          f"{function.rsplit('.', 1)[-1]}")
    return frames


def rtio_channel_names(device_db):
    """{RTIO channel: 'device/alias/...'} from a device_db. Satellite channels
    carry their destination in the upper bits there, as in RTIO errors."""
    names = {}
    for name, entry in device_db.items():
        if isinstance(entry, dict) and entry.get("type") == "local":
            channel = entry.get("arguments", {}).get("channel")
            if isinstance(channel, int):
                names.setdefault(channel, []).append(name)
    for alias, target in device_db.items():
        if not isinstance(target, str):
            continue
        seen = {alias}
        while isinstance(target, str) and target not in seen:
            seen.add(target)
            target = device_db.get(target)
        if isinstance(target, dict) and target.get("type") == "local":
            channel = target.get("arguments", {}).get("channel")
            if isinstance(channel, int):
                names.setdefault(channel, []).append(alias)
    return {channel: "/".join(found) for channel, found in names.items()}


class _GeneralVariableScanMasterSatelliteMixin(_DatasetRedirectMixin):
    """Shared implementation for the public master-satellite GVS files.

    This mixin deliberately does not inherit EnvExperiment, so ARTIQ Explorer
    cannot expose this common implementation as a third dashboard experiment.
    """

    VALID_NODES = ("Node1", "Node2")
    # The GUI argument of the same name overrides this per run; the class
    # value serves scans built without build() (repository examination, the
    # hand-built scans of the unit tests).
    pair_scan_sequences = False
    LEGACY_NODE_NAMES = {"Node1": "alice", "Node2": "bob"}
    EXPERIMENT_MODE = None
    # Overridden by tests to inject a stand-in for AOMPowerStabilizer; None
    # means Base builds the real one.
    _stabilizer_factory = None

    def _build_master_satellite_scan(self):
        if self.EXPERIMENT_MODE not in ("single_node", "two_nodes"):
            raise RuntimeError(
                "Public master-satellite GVS must define a fixed "
                "EXPERIMENT_MODE."
            )
        self.experiment_mode = self.EXPERIMENT_MODE
        # Repository examination intentionally supplies None for argument
        # values. Bind the suffixed two-node device superset now; real mode
        # validation, variable loading, and compatibility presentation happen
        # in prepare() after submitted values are available.
        self.base = BaseExperimentMasterSatellite(experiment=self)
        self.base.build()

        self.setattr_argument(
            "n_measurements",
            NumberValue(
                100,
                ndecimals=0,
                step=1,
                type="int",
            ),
        )
        self.setattr_argument(
            "scan_variable1_name",
            StringValue(
                "f_FORT" if self.EXPERIMENT_MODE == "single_node"
                else "f_FORT_Node1"
            ),
        )
        self.setattr_argument(
            "scan_sequence1",
            StringValue(
                "np.array([self.f_FORT])"
                if self.EXPERIMENT_MODE == "single_node"
                else "np.array([self.f_FORT_Node1])"
            ),
        )
        self.setattr_argument("scan_variable2_name", StringValue(""))
        self.setattr_argument(
            "scan_sequence2", StringValue("np.zeros(1)")
        )
        # Off: the full grid, every scan_sequence2 value at every
        # scan_sequence1 value. On: the two step together, point i using the
        # i-th value of each -- e.g. both nodes' target_780_HWP rotated
        # together in Two_nodes_atom_atom_parity_experiment.
        self.setattr_argument(
            "pair_scan_sequences",
            BooleanValue(False),
            tooltip="Step scan_variable2 together with scan_variable1 (point "
                    "i takes the i-th value of both sequences) instead of the "
                    "full grid. Needs scan_variable2_name and two equally "
                    "long sequences.",
        )
        # What both nodes share goes in the global dictionary, read in both
        # modes: the master-satellite globals and, in two_nodes mode, the
        # two-node sequence settings written bare ({'t_pumping': ...} =
        # t_pumping_two_nodes). A per-node name there raises
        # (merge_per_node_override_dictionaries). Then one dictionary per
        # node; both may stay populated at once, so switching selected_node
        # is the only edit needed to move a scan between nodes. A global in
        # a node's dictionary still works (saved submissions) but logs a
        # warning.
        self.setattr_argument(
            "override_ExperimentVariables_global",
            StringValue("{}"),
            tooltip="Overrides of variables both nodes share: the "
                    "master-satellite globals (two_atom_threshold, "
                    "t_delay_between_shots, ...) and, in two_nodes mode, the "
                    "two-node sequence settings written bare (t_pumping = "
                    "t_pumping_two_nodes, n_excitation_attempts, "
                    "t_OP_coil_relax_time, ...). Per-node variables are "
                    "refused here; they go in the node fields.",
        )
        for node in self.VALID_NODES:
            self.setattr_argument(
                f"override_ExperimentVariables_{node}",
                StringValue("{}"),
                tooltip=f"Overrides applied only when {node} runs. Otherwise "
                        f"ignored entirely -- not even evaluated, so it may "
                        f"name {node} variables that do not exist on the "
                        f"other node. Shared variables go in "
                        f"override_ExperimentVariables_global.",
            )
        # Per-point RTIOUnderflow retry, on by default since 2026-10-28, when
        # a retried point stopped keeping the failed attempt's data
        # (_run_scan_point). The three retry values below take effect only
        # when enable_Catch_UnderFlow is True.
        self.setattr_argument(
            "enable_Catch_UnderFlow",
            BooleanValue(True),
            "Catch Underflow",
        )
        self.setattr_argument(
            "underflow_max_retries",
            NumberValue(10, ndecimals=0, step=1, type="int"),
            "Catch Underflow",
        )
        self.setattr_argument(
            "underflow_backoff_ms",
            NumberValue(200.0, step=0.5),
            "Catch Underflow",
        )
        # Off by default: a skipped point leaves no measurements, so analysis
        # that maps measurements to scan points by position (n_measurements
        # per point) would label every later point wrongly. With it off the
        # run stops instead, and every point in the file is whole.
        self.setattr_argument(
            "skip_only_that_iteration_if_exhausted",
            BooleanValue(False),
            "Catch Underflow",
        )
        self.setattr_argument(
            "control_experiment",
            BooleanValue(False),
            "Control experiment",
        )

        self.experiment_function_registry = self._active_function_registry()
        function_names = sorted(self.experiment_function_registry)
        if not function_names:
            raise RuntimeError(
                f"No experiment functions are available for "
                f"{self.EXPERIMENT_MODE!r} mode."
            )
        self.setattr_argument(
            "experiment_function", EnumerationValue(function_names)
        )
        self.setattr_argument(
            "create_applets",
            BooleanValue(True),
            "Applets",
            tooltip="Ask the dashboard to show this node's applets, with "
                    "dataset names resolved for the selected node. Untick to "
                    "leave the dashboard exactly as it is. Requires the applet "
                    "dock's CCB policy to be 'Create and enable/disable "
                    "applets'.",
        )

    def _active_function_registry(self):
        if self.EXPERIMENT_MODE == "single_node":
            return build_single_node_function_registry()
        if self.EXPERIMENT_MODE == "two_nodes":
            return build_two_node_function_registry()
        raise RuntimeError("Invalid fixed master-satellite GVS mode.")

    def _publish_legacy_node_compatibility(self):
        if self.EXPERIMENT_MODE != "single_node":
            return
        try:
            self.which_node = self.LEGACY_NODE_NAMES[self.selected_node]
        except KeyError:
            raise ValueError(
                f"Unsupported selected_node {self.selected_node!r}; "
                "expected 'Node1' or 'Node2'."
            ) from None

    @staticmethod
    def _select_experiment_function(registry, name):
        try:
            return registry[name]
        except KeyError:
            available = ", ".join(sorted(registry)) or "<none>"
            raise ValueError(
                f"Experiment function {name!r} is not available in the "
                f"active execution mode. Available functions: {available}"
            ) from None

    @staticmethod
    def _evaluate_expression(expression, description, local_values):
        try:
            return eval(expression, globals(), local_values)
        except Exception as error:
            raise ValueError(
                f"Could not evaluate {description} {expression!r}: {error}"
            ) from error

    def _active_override_nodes(self):
        """Nodes whose per-node override dictionary applies to this run."""
        if self.EXPERIMENT_MODE == "single_node":
            return (self.selected_node,)
        return tuple(self.VALID_NODES)

    def _collect_authoritative_overrides(self):
        """Merge this run's global and per-node override dictionaries.

        The merge itself lives in the base module and is shared with
        MicrowaveScanOptimizer_master_satellite, so both experiments treat
        the idle node's dictionary identically -- see
        merge_per_node_override_dictionaries for why that dictionary is never
        even evaluated.

        globals() is passed so override text keeps being evaluated in THIS
        module's namespace, where entries like {'f_FORT': 240*MHz} find their
        units.
        """
        return merge_per_node_override_dictionaries(
            experiment=self,
            active_nodes=self._active_override_nodes(),
            resolve_target=self.base.resolve_experiment_variable_target,
            eval_globals=globals(),
            global_variable_names=self.base.global_variable_names,
        )

    def prepare(self):
        # n_measurements is a run-local GUI value sharing its name with a
        # persistent global dataset. Capture the submitted value before
        # configure_execution loads that dataset over the same attribute,
        # and re-assert it so the submitted value wins for this run.
        self.execution_n_measurements = self.n_measurements
        selected_node = (
            self.selected_node
            if self.EXPERIMENT_MODE == "single_node"
            else None
        )
        self.base.configure_execution(
            self.EXPERIMENT_MODE, which_node=selected_node
        )
        self.n_measurements = self.execution_n_measurements
        self._publish_legacy_node_compatibility()
        self.base.prepare()

        # The reused experiment_functions reach for the per-channel feedback
        # objects (self.stabilizer_FORT, self.stabilizer_AOM_A5, ...), which
        # only exist once AOMPowerStabilizer has been constructed -- it is what
        # publishes them onto the experiment. Without this, any function that
        # touches feedback fails to COMPILE, not merely to run. AOMsCoils and
        # the microwave optimizer already do this; two-node mode must not,
        # because master-satellite feedback is single-node only for now.
        # Both modes now. In two_nodes mode prepare_laser_stabilizer delegates
        # to prepare_laser_stabilizers_both_nodes, which builds one
        # AOMPowerStabilizer per crate writing the same per-node datasets a
        # single-node run writes -- so the applets and analysis notebooks are
        # unaffected. It used to be skipped here, which is why a two-node run
        # could only ever be open loop.
        self.base.prepare_laser_stabilizer(
            stabilizer_factory=self._stabilizer_factory
        )

        self.scan_variable1 = self.base.resolve_experiment_variable_target(
            str(self.scan_variable1_name)
        )
        self.scan_sequence1 = self._evaluate_expression(
            self.scan_sequence1,
            "scan_sequence1",
            {"self": self},
        )
        if len(self.scan_sequence1) == 0:
            raise ValueError("scan_sequence1 must contain at least one value.")

        scan_variable2_name = str(self.scan_variable2_name)
        if scan_variable2_name:
            self.scan_variable2 = (
                self.base.resolve_experiment_variable_target(
                    scan_variable2_name
                )
            )
            self.scan_sequence2 = self._evaluate_expression(
                self.scan_sequence2,
                "scan_sequence2",
                {"self": self},
            )
            if len(self.scan_sequence2) == 0:
                raise ValueError(
                    "scan_sequence2 must contain at least one value."
                )
        else:
            self.scan_variable2 = None
            self.scan_sequence2 = np.zeros(1)
        self._validate_paired_scan()

        self.authoritative_overrides = self._collect_authoritative_overrides()

        self.experiment_name = str(self.experiment_function)
        # Filename suffix built the same way standalone GeneralVariableScan
        # does it (GeneralVariableScan.py:137-140, scan vars joined with
        # "_and_"), so results from both stacks stay filterable by the same
        # substrings in the Analysis notebooks.
        scan_vars = [
            name for name in (str(self.scan_variable1_name),
                              str(self.scan_variable2_name)) if name
        ]
        self.scan_var_filesuffix = "_and_".join(scan_vars)
        self._selected_experiment_function = (
            self._select_experiment_function(
                self.experiment_function_registry,
                self.experiment_name,
            )
        )
        self._initialize_legacy_result_attributes()
        self.needs_experiment_variable_reload = (
            self._has_earlier_queued_experiment()
        )

    def _initialize_legacy_result_attributes(self):
        self.measurement = 0
        for detector in (
            "SPCM0",
            "SPCM1",
            "SPCM0_OtherNode",
            "SPCM1_OtherNode",
        ):
            setattr(self, f"{detector}_RO1", 0)
            setattr(self, f"{detector}_RO2", 0)

    def _has_earlier_queued_experiment(self):
        status = self.scheduler.get_status()
        rid = self.scheduler.rid
        earlier_count = len(
            [scheduled_rid for scheduled_rid in status if scheduled_rid < rid]
        )
        logging.info(
            "RID %s has %s earlier queued experiment(s)", rid, earlier_count
        )
        return earlier_count > 0

    def _refresh_runtime_variable_state(self):
        self.base.refresh_compatibility_variables()
        self.base.refresh_variable_dependent_state()

    def _apply_run_wide_overrides(self):
        # A scanned variable takes its value from the scan at every point, so
        # an override of it never actually runs; recording one would put a
        # value in the h5 that was never used.
        scanned = {getattr(self, "scan_variable1", None),
                   getattr(self, "scan_variable2", None)}
        for target, value in self.authoritative_overrides.items():
            setattr(self, target, value)
            if target in scanned:
                continue
            # Record the value this run actually used, as the standalone
            # GeneralVariableScan has always done (initialize_datasets:
            # self.set_dataset(var, val)). The h5 then holds the persistent
            # value under `archive` and this one under `datasets`, and
            # Analysis loads datasets after archive, so a notebook sees what
            # ran. Without it, RID 38903 -- a Z-parity run, t_MW_RF_pulse
            # overridden to 0 -- was analysed as X parity from the archived
            # 69.6 us.
            #
            # LOCAL to this run, and spelled out: broadcast or persist would
            # write the master's dataset DB, making the override the new
            # default and clearing the persist flag -- the trap that once
            # de-persisted n_measurements.
            #
            # `target` is the name the run read and archived: suffixed for a
            # per-node variable (t_MW_RF_pulse_Node2), bare for a global. The
            # bare per-node name is the STANDALONE stack's dataset, which
            # this run never read, so it is not written.
            #
            # ORDER MATTERS. A local write to a key this run has already
            # BROADCAST makes ARTIQ delete that key from the master DB
            # (artiq/master/worker_db.py, DatasetManager.set: the
            # `elif key in self._broadcaster.raw_view` branch). This pass runs
            # before initialize_result_state, which broadcasts and persists
            # n_measurements, and before any other broadcast of a variable
            # name -- verified by replaying RID 38903 through the real
            # DatasetManager, 2026-10-14 -- so nothing can be deleted. Keep it
            # first in _initialize_run_state;
            # test_override_records_precede_every_broadcast pins that.
            self.set_dataset(target, value, broadcast=False, persist=False)
        self._refresh_runtime_variable_state()

    def _apply_scan_point(self, variable1_value, variable2_value=None):
        setattr(self, self.scan_variable1, variable1_value)
        if self.scan_variable2 is not None:
            setattr(self, self.scan_variable2, variable2_value)
        self._refresh_runtime_variable_state()

    @kernel
    def initialize_hardware(self):
        self.base.initialize_hardware()

    @kernel
    def force_other_node_off(self):
        self.base.force_other_node_off()

    def _initialize_run_state(self):
        """Refresh queued values and initialize transient run results once."""
        if self.needs_experiment_variable_reload:
            self.base.reload_experiment_variables()
            # n_measurements is a run-local GUI value, not a persistent write.
            self.n_measurements = self.execution_n_measurements

        self._apply_run_wide_overrides()
        self.base.initialize_result_datasets()

        # initialize_result_datasets() only covers progress, the
        # _current_iteration histograms and the magnetometer series.
        # The atom-physics functions also need the full result surface (SPCM
        # datasets, per-measurement buffers and the host scalars they read in
        # kernels), which the microwave optimizer already sets up this way.
        #
        # Both modes now. The names are the same unsuffixed ones either way,
        # because a two-node run's results describe the JOINT measurement --
        # one atom per trap, read out through one gate of the four
        # master-local SPCMs -- so there is nothing per-node to suffix. Until
        # 2026-10-06 this was single_node-only and the method refused outside
        # it, which is what stopped any two-node sequence from compiling: the
        # host scalars have to EXIST before compilation, not merely before
        # running.
        self.base.initialize_result_state()

        # Scan labels, exactly as the standalone GeneralVariableScan publishes
        # them; the retention applet reads these for its axis.
        #
        # Both modes. These used to be single_node-only on the grounds that
        # Base publishes scan_var_dataset as part of the legacy-compatibility
        # namespace -- but it is now published unsuffixed in both modes,
        # because it names run-global state rather than a node's wiring. The
        # omission mattered more than it looks: a scan point mutates
        # ATTRIBUTES, not datasets, so the archived per-node values are the
        # pre-scan ones and these three were the only record of what a
        # two-node run actually scanned. There was none.
        scan_names = [
            name for name in (str(self.scan_variable1_name),
                              str(self.scan_variable2_name)) if name
        ]
        self.set_dataset(self.scan_var_dataset, ",".join(scan_names),
                         broadcast=True)
        self.set_dataset(self.scan_sequence1_dataset, self.scan_sequence1,
                         broadcast=True)
        self.set_dataset(self.scan_sequence2_dataset, self.scan_sequence2,
                         broadcast=True)

        self._create_node_applets()

    def _create_node_applets(self):
        """Ask the dashboard for applets, if enabled.

        Two-node mode gets the SHARED applets only. create_applets_for
        refuses anything but single_node, because the per-node applets read
        result datasets that two-node mode does not suffix, so both nodes
        would write the same names and there would be nothing for a per-node
        applet to point at. The shared specs do not have that problem: they
        name their datasets explicitly per node, so they resolve unchanged in
        either mode.

        Feedback RF works in two-node mode too: each node's stabilizer
        writes p_AOM_A1_history_NodeX, the name a single-node run writes, so
        the shared feedback RF applets plot live values in either mode.

        Never fatal: the applets are a convenience, and the CCB is a
        dashboard-side service that artiq_run and the offline compile check
        do not provide meaningfully.
        """
        if not self.create_applets:
            return
        try:
            if self.EXPERIMENT_MODE == "single_node":
                from applets_master_satellite import (
                    APPLET_SPECS,
                    GVS_APPLET_SPECS,
                    SCAN_APPLET_SPECS,
                    SHARED_APPLET_SPECS,
                    create_applets_for,
                )

                # GVS_APPLET_SPECS is the atom-check loading histogram, asked
                # for here rather than living in APPLET_SPECS: the dataset is
                # seeded on every experiment, so a default entry would open an
                # empty histogram on runs that never load an atom.
                #
                # SCAN_APPLET_SPECS ("Time without atom") goes through
                # shared_specs instead, so GVS and MicrowaveScanOptimizer both
                # create the SAME applet rather than one each under different
                # groups. See its note in applets_master_satellite.
                created = create_applets_for(
                    self, self.base,
                    specs=APPLET_SPECS + GVS_APPLET_SPECS,
                    shared_specs=SHARED_APPLET_SPECS + SCAN_APPLET_SPECS,
                )
            else:
                from applets_master_satellite import (
                    FIRST_SHOT_ALIVE_APPLET_SPECS,
                    FIRST_SHOT_ALIVE_EXPERIMENTS,
                    SCAN_APPLET_SPECS,
                    SHARED_APPLET_SPECS,
                    TWO_NODE_ALTERNATING_APPLET_SPECS,
                    TWO_NODE_ALTERNATING_EXPERIMENTS,
                    TWO_NODE_APPLET_SPECS,
                    TWO_NODE_PER_NODE_RETENTION_APPLET_SPECS,
                    TWO_NODE_PER_NODE_RETENTION_EXPERIMENTS,
                    TWO_NODE_PER_NODE_RETENTION_RETIRES,
                    create_two_node_applets_for,
                )

                # The per-node atom-check histograms are asked for only by the
                # experiments whose loader actually writes those datasets.
                # Putting them in TWO_NODE_APPLET_SPECS would open two empty
                # histograms on every run that uses the joint loader instead --
                # the same reason GVS_APPLET_SPECS exists rather than growing
                # APPLET_SPECS.
                specs = TWO_NODE_APPLET_SPECS
                retire_titles = ()
                if self.experiment_name in TWO_NODE_ALTERNATING_EXPERIMENTS:
                    specs = specs + TWO_NODE_ALTERNATING_APPLET_SPECS
                # The node-resolved SHOTS REPLACE the joint plots rather than
                # adding to them: this experiment never calls first_shot or
                # second_shot, so AllSPCMs_RO1/RO2 stay at zero all run. So
                # swap in the six per-node plots and close the three joint
                # ones -- creating the joint retention applet here only to
                # disable it a line later would be worse than not asking.
                if self.experiment_name in TWO_NODE_PER_NODE_RETENTION_EXPERIMENTS:
                    specs = TWO_NODE_PER_NODE_RETENTION_APPLET_SPECS
                    retire_titles = TWO_NODE_PER_NODE_RETENTION_RETIRES
                # The time each kept measurement really took to load when the
                # first shot has to confirm the atoms. Only the experiments
                # whose confirm loop writes it, and only with the flag on;
                # otherwise it is not created, and so it is closed.
                if (self.advance_if_atom_in_first_shot
                        and self.experiment_name in FIRST_SHOT_ALIVE_EXPERIMENTS):
                    specs = specs + FIRST_SHOT_ALIVE_APPLET_SPECS

                # NOT create_shared_applets_for: that deliberately runs no
                # retirement pass, so a single-node run's "retention and
                # loading AllSPCMs" would stay open through this run and plot
                # the same live AllSPCMs_RO1/RO2 against single_atom_threshold.
                created = create_two_node_applets_for(
                    self, self.base,
                    specs=specs,
                    shared_specs=SHARED_APPLET_SPECS + SCAN_APPLET_SPECS,
                    retire_titles=retire_titles,
                )
        except Exception as error:  # noqa: BLE001 - convenience only
            logging.warning("could not create applets: %s", error)
        else:
            # which_node is None in two-node mode; name the mode instead.
            logging.info("requested %d applets for %s", len(created),
                         self.base.which_node or self.EXPERIMENT_MODE)

    def _execute_scan_point(self, variable1_value, variable2_value, iteration):
        """Execute one scan point without rebuilding or preparing devices."""
        self._apply_scan_point(variable1_value, variable2_value)
        self.initialize_hardware()
        self.base.reset_result_state_for_scan_point()
        self.set_dataset("iteration", iteration, broadcast=True)
        self._selected_experiment_function(self)
        self._save_scan_point_results()

    def _save_scan_point_results(self):
        # Saved per scan point, matching standalone GeneralVariableScan's
        # run_iteration: it lets a run be quit early without losing data, and
        # guards against ARTIQ corrupting the h5 during worker cleanup (see
        # README on write_results). The master-satellite port had no named
        # save at all, so GVS results carried no node at all -- ARTIQ's own
        # filename is just <rid>-GeneralVariableScan_master_satellite_
        # single_node.h5. result_name_tag() puts the node first in the name.
        self.write_results({
            'name': self.base.result_name_tag() + "_"
                    + self.experiment_name[:-11]
                    + "_scan_over_" + self.scan_var_filesuffix
        })

    def _validate_paired_scan(self):
        """A paired scan needs a second variable and equally long sequences."""
        if not self.pair_scan_sequences:
            return
        if self.scan_variable2 is None:
            raise ValueError(
                "pair_scan_sequences steps scan_variable2 with "
                "scan_variable1: set scan_variable2_name too.")
        if len(self.scan_sequence2) != len(self.scan_sequence1):
            raise ValueError(
                f"pair_scan_sequences needs equally long sequences: "
                f"scan_sequence1 has {len(self.scan_sequence1)} values, "
                f"scan_sequence2 {len(self.scan_sequence2)}.")

    def _scan_points(self):
        """(variable1, variable2) per scan point, in run order."""
        if self.pair_scan_sequences:
            return list(zip(self.scan_sequence1, self.scan_sequence2))
        return [(variable1_value, variable2_value)
                for variable1_value in self.scan_sequence1
                for variable2_value in self.scan_sequence2]

    def _report_underflow(self, message):
        logging.warning(message)
        print(message)

    def _underflow_context(self, snapshot):
        """What a failed attempt got through, read from its own appends before
        roll_back_list_datasets cuts them: measurements finished, measurements
        past the first-shot confirm (one more than finished = it failed inside
        a measurement's excitation loop or second shot), and its last
        atom_loading_wall_clock stamp in mu."""
        def appended(name):
            entry = snapshot.get(self.base.resolve_result_dataset_name(name))
            if entry is None:
                return []
            value, length = entry
            return value[length:]

        loads = appended("atom_loading_wall_clock")
        return {
            "measurements": len(appended("AllSPCMs_RO1")),
            "confirmed": len(appended("n_first_shot_confirm_tries")),
            "last_load_mu": int(loads[-1]) if loads else None,
        }

    def _rtio_channel_name(self, channel):
        if self._rtio_channel_names is None:
            try:
                self._rtio_channel_names = rtio_channel_names(
                    self.get_device_db())
            except (AttributeError, KeyError, TypeError):
                self._rtio_channel_names = {}
        return self._rtio_channel_names.get(channel, "unknown device")

    def _record_underflow(self, error, iteration, retry, context):
        """Keep each underflow in the h5: the rollback leaves nothing of the
        failed attempt and the warning reaches only the dashboard log.

        underflow_channel, _timestamp_mu (-1: message not parsed),
        _slack_mu (NaN for a DMA playback), _measurements_done,
        _confirmed_measurements, _mu_after_last_load (NaN: no load in the
        attempt) per underflow; underflow_sites (one line each) and
        underflow_messages (message and core device traceback) as text.

        Returns (one-line summary, core device traceback or '').
        """
        parsed = parse_underflow_message(str(error))
        timestamp, channel, slack = parsed if parsed else (-1, -1, None)
        core_exception = getattr(error, "artiq_core_exception", None)
        traceback = str(core_exception) if core_exception is not None else ""
        frames = repository_frames(traceback)
        last_load = context["last_load_mu"]
        after_load = (timestamp - last_load
                      if parsed and last_load is not None else None)

        if parsed is None:
            kind = "unparsed message"
        else:
            kind = (("DMA playback" if slack is None else f"CPU write, slack {slack} mu")
                    + f", channel {channel:#x} = {self._rtio_channel_name(channel)}")
        summary = (
            f"{kind}; at {', '.join(frames) if frames else 'unknown line'}; "
            f"the attempt had finished {context['measurements']} measurement(s), "
            f"{context['confirmed']} past the first-shot confirm"
            + ("" if after_load is None
               else f", {after_load / 1e6:.3f} ms after its last load stamp"))

        events = self._underflow_events
        events["underflow_iteration"].append(iteration)
        events["underflow_retry"].append(retry)
        events["underflow_channel"].append(channel)
        events["underflow_timestamp_mu"].append(timestamp)
        events["underflow_slack_mu"].append(
            float("nan") if slack is None else float(slack))
        events["underflow_measurements_done"].append(context["measurements"])
        events["underflow_confirmed_measurements"].append(context["confirmed"])
        events["underflow_mu_after_last_load"].append(
            float("nan") if after_load is None else float(after_load))
        self._underflow_sites.append(
            f"iteration {iteration} retry {retry}: {summary}")
        self._underflow_messages.append(
            f"iteration {iteration} retry {retry}: {error}\n{traceback}")
        # fresh lists every time, so a later attempt's rollback leaves them be
        for key, values in events.items():
            self.set_dataset(key, list(values), broadcast=False, persist=False)
        self.set_dataset("underflow_sites", "\n".join(self._underflow_sites),
                         broadcast=False, persist=False)
        self.set_dataset("underflow_messages",
                         "\n\n".join(self._underflow_messages),
                         broadcast=False, persist=False)
        return summary, traceback

    def _underflow_backoff(self):
        """Pause on the host; the retried point resets the core via Base."""
        time.sleep(max(0.0, float(self.underflow_backoff_ms)) / 1000.0)

    def _record_underflow_retries(self, iteration, retries):
        """underflow_retries_per_scan_point: the failed attempts of each point.

        Set before every attempt, so the attempt that succeeds saves its
        point's final count with the point's results.
        """
        self._underflow_retries_per_scan_point[iteration:] = [retries]
        self.set_dataset(
            "underflow_retries_per_scan_point",
            list(self._underflow_retries_per_scan_point),
            broadcast=False,
            persist=False,
        )

    def _run_scan_point(self, variable1_value, variable2_value, iteration):
        """Run one scan point, retrying only on RTIOUnderflow when enabled.

        With enable_Catch_UnderFlow off the point executes directly and any
        failure propagates unchanged. When it is on (the default), only
        RTIOUnderflow is caught; a retried point repeats the full scan-point
        flow, including hardware initialization and its core reset. DRTIO,
        SPI, resolution, and all other failures always propagate.

        A failed attempt leaves no data: before the retry, every list dataset
        is cut back to its length before the attempt
        (roll_back_list_datasets). Until 2026-10-28 the attempt's
        measurements stayed in the series ahead of the retry's, so anything
        that maps measurements to scan points by position (n_measurements per
        point, as the analysis notebooks do) read every later point with the
        wrong scan value -- RID 39594 held 502 measurements for 5 x 100. Now
        each point holds exactly its one successful attempt, or nothing when
        the retries run out and skip_only_that_iteration_if_exhausted is on.
        The cut happens before the run stops, too, so the final h5 has no
        partial point either.

        Recorded: underflow_retries_per_scan_point (failed attempts per
        point) and underflow_skipped_scan_points (iterations given up on).
        """
        if not self.enable_Catch_UnderFlow:
            self._execute_scan_point(
                variable1_value, variable2_value, iteration
            )
            return

        dataset_mgr = self._HasEnvironment__dataset_mgr
        retries = 0
        while True:
            self._record_underflow_retries(iteration, retries)
            snapshot = snapshot_list_datasets(dataset_mgr)
            try:
                self._execute_scan_point(
                    variable1_value, variable2_value, iteration
                )
                return
            except RTIOUnderflow as error:
                ### read what the attempt got through before the cut
                context = self._underflow_context(snapshot)
                removed = roll_back_list_datasets(dataset_mgr, snapshot)
                retries += 1
                self._record_underflow_retries(iteration, retries)
                maximum_retries = int(self.underflow_max_retries)
                summary, traceback = self._record_underflow(
                    error, iteration, retries, context)
                self._report_underflow(
                    f"RTIO underflow at iteration {iteration}, "
                    f"retry {retries}/{maximum_retries}: {error}. {summary}. "
                    f"Discarded the attempt's data: "
                    f"{_describe_rollback(removed)}"
                )
                if traceback:
                    self._report_underflow(traceback)
                if retries >= maximum_retries:
                    message = (
                        f"RTIO underflow at iteration {iteration} exceeded "
                        f"max retries ({maximum_retries})."
                    )
                    self._report_underflow(message)
                    if not self.skip_only_that_iteration_if_exhausted:
                        raise
                    self._underflow_skipped_scan_points.append(iteration)
                    self.set_dataset(
                        "underflow_skipped_scan_points",
                        list(self._underflow_skipped_scan_points),
                        broadcast=False,
                        persist=False,
                    )
                    self._save_scan_point_results()
                    return
                self._underflow_backoff()

    def run(self):
        self._initialize_run_state()
        self._underflow_retries_per_scan_point = []
        self._underflow_skipped_scan_points = []
        self._underflow_events = {key: [] for key in (
            "underflow_iteration", "underflow_retry", "underflow_channel",
            "underflow_timestamp_mu", "underflow_slack_mu",
            "underflow_measurements_done", "underflow_confirmed_measurements",
            "underflow_mu_after_last_load")}
        self._underflow_sites = []
        self._underflow_messages = []
        self._rtio_channel_names = None
        if self.enable_Catch_UnderFlow:
            # present, if empty, in every file of a run with the retry on
            self.set_dataset("underflow_skipped_scan_points", [],
                             broadcast=False, persist=False)

        # Once per run, before the scan loop. single_node mode leaves the
        # other node out of the initialization lifecycle, so without this its
        # beams and coils keep whatever state the previous run left them in.
        # Placing it here is safe: the core.reset() inside the per-scan-point
        # initialize_hardware() can only drive TTLs low and cannot change a
        # held Zotino output, so it cannot undo this. The host-side mode guard
        # keeps the kernel out of the two-node variant, which has no "other"
        # node and whose device lists would not type-check.
        if self.EXPERIMENT_MODE == "single_node":
            self.force_other_node_off()

        iteration = 0
        self.set_dataset("iteration", iteration, broadcast=True)
        for variable1_value, variable2_value in self._scan_points():
            self._run_scan_point(variable1_value, variable2_value, iteration)
            iteration += 1

        print(
            "**************** General Variable Scan master-satellite DONE "
            "****************"
        )
