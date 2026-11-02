"""Shared scan implementation (mixins) for the master-satellite GVS files.

This module holds no public EnvExperiment class and is not an ARTIQ Explorer
experiment; the runnable experiments live in the mode-named GVS files.
"""

import copy
import inspect
import logging
import re
import time
import types

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


# ---------------------------------------------------------------------------
# Microwave health check (run_health_check; since 2026-10-29)
#
# The FORT polarization drifts, so the AC Stark shift and with it the
# microwave transition frequencies drift during a long scan. Every
# health_check_every_n_ite scan points the GVS hands each selected transition
# to MicrowaveScanOptimizer_master_satellite -- the experiment that already
# knows how to check and re-centre one -- and WAITS for it inside its own run,
# so the scan stays one run and one h5 file. See _run_health_check_event.
# ---------------------------------------------------------------------------

HEALTH_CHECK_OPTIMIZER_FILE = (
    "qn_artiq_routines\\MicrowaveScanOptimizer_master_satellite.py")
HEALTH_CHECK_OPTIMIZER_CLASS = "MicrowaveScanOptimizer_master_satellite"
#: The optimizer resubmits its own scan at 99 after a failed check, so the
#: whole chain runs at 99; a GVS at 99 itself would never yield to it.
HEALTH_CHECK_PRIORITY = 99

#: label, the optimizer's scan flag (also the GVS argument name), and the
#: center that scan checks and re-writes -- scan_dict[flag]["center"] of the
#: optimizer, pinned by a test there (the optimizer module is not imported
#: here: it pulls scipy and the fit models into every GVS worker). IN CHECK
#: ORDER: 11 and m11 pulse 01 first and m10 pulses 00 first, so 00 and 01 are
#: re-centred before the transitions that depend on them.
HEALTH_CHECK_SCANS = (
    ("00", "Frequency_00_Scan", "f_microwaves_00_dds"),
    ("01", "Frequency_01_Scan", "f_microwaves_01_dds"),
    ("11", "Frequency_11_Scan", "f_microwaves_11_dds"),
    ("m10", "Frequency_m10_Scan", "f_microwaves_m10_dds"),
    ("m11", "Frequency_m11_Scan", "f_microwaves_m11_dds"),
)
HEALTH_CHECK_LABELS = tuple(label for label, _, _ in HEALTH_CHECK_SCANS)

#: Every scan flag the optimizer declares (its scan_options, pinned by a test
#: there). All of them go in the expid, False but the one checked, so a
#: changed optimizer default cannot add a second scan.
HEALTH_CHECK_OPTIMIZER_SCAN_FLAGS = (
    "Frequency_00_Scan", "Frequency_01_Scan", "Frequency_11_Scan",
    "Frequency_m10_Scan", "Frequency_m11_Scan",
    "Time_00_Scan", "Time_01_Scan", "Time_11_Scan", "Time_m10_Scan",
    "Time_m11_Scan",
    "Ramsey_00_Scan", "Ramsey_01_Scan", "Ramsey_11_Scan",
)

HEALTH_CHECK_NODE_SELECTIONS = {
    "Node1": ("Node1",),
    "Node2": ("Node2",),
    "both_nodes": ("Node1", "Node2"),
}

# ARTIQ 7 RunStatus names (artiq/master/scheduler.py). Not started = still
# queued (a Terminate leaves it there); finished = no longer to be waited for.
_HEALTH_CHECK_NOT_STARTED = frozenset(
    {"pending", "flushing", "preparing", "prepare_done"})
_HEALTH_CHECK_FINISHED = frozenset({"run_done", "analyzing", "deleting"})

# The calibrations the optimizer measures or writes. Never forwarded: the
# chain must start from, and write, the stored values.
_MICROWAVE_CALIBRATION = re.compile(
    r"^(f_microwaves_\w*dds|t_microwave_\w*pulse|t_MW_RF_pulse|health_check_uw_)")


def republish_run_datasets(dataset_mgr, persistent_keys=("n_measurements",)):
    """Publish this run's broadcast datasets to the master again, unchanged.

    For after another experiment ran inside this one (the health check): the
    optimizer broadcasts the same result names -- AllSPCMs_*, the scan labels,
    the histograms, n_measurements -- so the master DB and every applet show
    ITS values until this run writes each name again. This run's own copies
    (and its h5) are untouched by that; only the master's view is.

    Every NON-persistent broadcast key is set again, plus *persistent_keys*:
    n_measurements is persisted by every experiment, and the last optimizer
    run's n_measurements (the check's or its scan's) would otherwise be the
    stored default and the retention applet's grouping. No other persistent key: feedback results
    (p_*_NodeX) and calibrations are the latest MEASURED values, and the
    optimizer's are newer than this run's.

    Same object, same flags, the way roll_back_list_datasets does it: a
    broadcast key must stay broadcast, since a local write deletes it from the
    master DB. Returns the keys published.
    """
    broadcast_view = dataset_mgr._broadcaster.raw_view
    published = []
    for key, (persist, value) in list(broadcast_view.items()):
        if persist and key not in persistent_keys:
            continue
        dataset_mgr.set(key, value, broadcast=True, persist=persist,
                        archive=key in dataset_mgr.local)
        published.append(key)
    return published


def _override_literal(value):
    """Python text for one forwarded override value, evaluated by the
    optimizer in its own module namespace (numpy's array is imported there).
    numpy scalars become Python numbers, and a non-finite float is spelled
    float('nan'), which eval understands and the bare repr 'nan' is not."""
    if isinstance(value, np.ndarray):
        return f"array({_override_literal(value.tolist())})"
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return f"float({str(value)!r})"
    if isinstance(value, list):
        return "[" + ", ".join(_override_literal(v) for v in value) + "]"
    if isinstance(value, tuple):
        inner = ", ".join(_override_literal(v) for v in value)
        return "(" + inner + ("," if len(value) == 1 else "") + ")"
    return repr(value)


def format_override_dictionary(overrides):
    """The text of an override_ExperimentVariables_<node> argument."""
    return "{" + ", ".join(
        f"{name!r}: {_override_literal(value)}"
        for name, value in overrides.items()) + "}"


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
    # Likewise for the microwave health check: off unless submitted on.
    run_health_check = False
    # Seconds between scheduler polls while a health check waits for a chain
    # that is not ready to run yet; tests set 0.
    health_check_poll_interval = 0.5
    # Seconds to wait before resubmitting a chain that ended with no result
    # of its own (the optimizer crashed); tests set 0.
    health_check_crash_backoff_s = 60
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
                    "leave the dashboard exactly as it is -- except the 'MW "
                    "health check status' applet, which run_health_check "
                    "always requests. Requires the applet "
                    "dock's CCB policy to be 'Create and enable/disable "
                    "applets'.",
        )
        self._build_health_check_arguments()

    def _build_health_check_arguments(self):
        """The "Health Check" group: all inert unless run_health_check.

        Names chosen to collide with no ExperimentVariable (a collision would
        need the capture-and-re-assert dance n_measurements needs).
        """
        group = "Health Check"
        self.setattr_argument(
            "run_health_check",
            BooleanValue(False),
            group,
            tooltip="Every health_check_every_n_ite scan points (never before "
                    "the first), pause and run MicrowaveScanOptimizer_master_"
                    "satellite on each ticked transition of the selected "
                    "node(s), one at a time, re-optimizing until it passes. "
                    "Needs priority below 99 and the main pipeline.",
        )
        self.setattr_argument(
            "health_check_every_n_ite",
            NumberValue(1, ndecimals=0, step=1, type="int"),
            group,
        )
        # A default although the choice is the user's: ARTIQ raises
        # DefaultMissing in build() for any expid without the key, which would
        # break every saved or scripted (artiq_client) GVS submission. The
        # dashboard shows the first choice either way.
        self.setattr_argument(
            "health_check_node_selection",
            EnumerationValue(list(HEALTH_CHECK_NODE_SELECTIONS),
                             default="Node1"),
            group,
            tooltip="Either node, or both, in either GVS mode.",
        )
        for _, flag, _ in HEALTH_CHECK_SCANS:
            self.setattr_argument(flag, BooleanValue(False), group)
        self.setattr_argument(
            "target_fidelity",
            NumberValue(0.80, ndecimals=2, step=0.01),
            group,
            tooltip="One target for every check.",
        )
        self.setattr_argument(
            "health_check_n_measurements",
            NumberValue(200, ndecimals=0, step=1, type="int"),
            group,
            tooltip="The check run's n_measurements.",
        )
        self.setattr_argument(
            "health_check_optimization_n_measurements",
            NumberValue(100, ndecimals=0, step=1, type="int"),
            group,
            tooltip="n_measurements of each point of the re-optimization "
                    "scan a failed check submits (the optimizer's "
                    "n_measurements_optimization_scan; 0 = the check's).",
        )
        self.setattr_argument(
            "health_check_adaptive_half_range_linewidths",
            NumberValue(1.0, ndecimals=2, step=0.5),
            group,
            tooltip="The optimizer's adaptive_half_range_linewidths.",
        )
        self.setattr_argument(
            "health_check_adaptive_points_per_linewidth",
            NumberValue(4.0, ndecimals=1, step=1),
            group,
            tooltip="The optimizer's adaptive_points_per_linewidth.",
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
        self._prepare_health_check()

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

    # ------------------------------------------------------------------
    # Microwave health check -- "pause and delegate" (user, 2026-10-29)
    #
    # Before scan point i (i > 0, i % health_check_every_n_ite == 0) and
    # OUTSIDE _run_scan_point, so an underflow retry never repeats it: for
    # each selected node, then each ticked transition in HEALTH_CHECK_SCANS
    # order, one at a time --
    #   1. submit MicrowaveScanOptimizer_master_satellite at 99 in this
    #      pipeline, parent_rid = this RID (_health_check_expid);
    #   2. release the core and pause() until every run of that chain (the
    #      check and the scan it submits when the check fails) is done;
    #   3. passed? A failed check: submit the SAME expid again at once, at 99
    #      so it runs next, without limit -- the growing count on the status
    #      applet is the signal to relock lasers. A chain with no result of
    #      its own (the optimizer crashed): the same, but only after
    #      health_check_crash_backoff_s;
    #   4. re-read that one center from the master DB.
    # Then put this run's datasets back on the dashboard and save. Terminate
    # ends the scan cleanly and deletes nothing from the queue.
    # ------------------------------------------------------------------

    def _health_check_selected_scans(self):
        return tuple(scan for scan in HEALTH_CHECK_SCANS
                     if getattr(self, scan[1], False))

    def _health_check_node_is_active(self, node):
        """Whether this run loaded, and so uses, the node's variables."""
        return (self.EXPERIMENT_MODE == "two_nodes"
                or node == getattr(self, "selected_node", None))

    def _prepare_health_check(self):
        """Refuse, in prepare, every setting the health check cannot work with.

          * Priority 99 or above. pause() yields only to a prepared run with
            a higher priority key, so a GVS at 99 never lets its own 99 chain
            run and the wait would spin for ever.
          * A pipeline other than main: the chain goes to main, and a run
            pauses only for runs of its own pipeline.
          * A selected center this run overrides or scans. The run's own copy
            would hide the optimizer's new value -- the point of the check.
          * Overrides the optimizer would refuse (_validate_forwarded_overrides);
            otherwise every check would crash in its prepare and be
            resubmitted for ever.
          * health_check_n_measurements below 1: every check would crash in
            its run (it divides by n_measurements) and be resubmitted for
            ever, every health_check_crash_backoff_s.
        """
        self._health_check_every_n = 0
        self._health_check_nodes = ()
        self._health_check_scans = ()
        if not self.run_health_check:
            return

        every_n = int(self.health_check_every_n_ite)
        if every_n < 1:
            raise ValueError(
                f"health_check_every_n_ite must be at least 1, got {every_n}.")
        check_n = int(self.health_check_n_measurements)
        if check_n < 1:
            raise ValueError(
                f"health_check_n_measurements must be at least 1, got "
                f"{check_n}.")
        priority = int(self.scheduler.priority)
        if priority >= HEALTH_CHECK_PRIORITY:
            raise ValueError(
                f"run_health_check needs a priority below "
                f"{HEALTH_CHECK_PRIORITY}, got {priority}: the optimizer chain "
                f"runs at {HEALTH_CHECK_PRIORITY}, and a run pauses only for a "
                f"higher priority, so this scan would wait for it for ever.")
        pipeline = self.scheduler.pipeline_name
        if pipeline != "main":
            raise ValueError(
                f"run_health_check needs the 'main' pipeline, got "
                f"{pipeline!r}: the optimizer chain is submitted to main, and "
                f"a run pauses only for runs of its own pipeline.")
        scans = self._health_check_selected_scans()
        if not scans:
            raise ValueError(
                "run_health_check is on but no transition is ticked: tick at "
                "least one of " + ", ".join(
                    flag for _, flag, _ in HEALTH_CHECK_SCANS) + ".")
        selection = str(self.health_check_node_selection)
        try:
            nodes = HEALTH_CHECK_NODE_SELECTIONS[selection]
        except KeyError:
            raise ValueError(
                f"Unknown health_check_node_selection {selection!r}; expected "
                f"one of {', '.join(HEALTH_CHECK_NODE_SELECTIONS)}.") from None

        hidden = set(self.authoritative_overrides) | {
            self.scan_variable1, self.scan_variable2}
        for node in nodes:
            for _, flag, center in scans:
                name = f"{center}_{node}"
                if name in hidden:
                    raise ValueError(
                        f"run_health_check with {flag} on {node}: this run "
                        f"overrides or scans {name}, so it would keep its own "
                        f"value and never use the one the check writes. "
                        f"Remove it from the overrides/scan, or untick "
                        f"{flag}.")
        for node in nodes:
            self._validate_forwarded_overrides(node)

        self._health_check_every_n = every_n
        self._health_check_nodes = nodes
        self._health_check_scans = scans

    def _health_check_forwarded_overrides(self, node):
        """What the optimizer on *node* should run with: {name: value}.

        The user's rule: "the shared global values should be overwritten in
        the optimizer as well because they depend on that". The optimizer
        always runs single_node and has no global field, but its node field
        takes bare names (resolved to that node), that node's suffixed names
        and plain globals. So:

          * two_nodes: ALL 13 shared two-node settings, written BARE with
            this run's CURRENT <bare>_two_nodes value -- stored or overridden
            alike, so the optimizer's per-node twin (t_pumping_Node1, ...)
            takes the value the two-node sequence runs with. A SCANNED one
            gets the run's value from before the scan
            (_health_check_shared_values), not the current point's;
          * this run's overrides of the node's own variables (two_nodes; in
            single_node only when it is the selected node -- there are no
            others). A two-node run cannot override a shared setting's
            per-node twin (t_pumping_Node1): its resolver refuses it, since
            the two-node sequence reads only <bare>_two_nodes. So the bare
            entries above never meet a twin entry ("overridden twice");
          * this run's overrides of plain globals.

        Never n_measurements (the optimizer has its own), a scan variable
        (its value is the current point's, not a setting), or a microwave
        calibration the optimizer measures or writes (_MICROWAVE_CALIBRATION):
        the chain must start from, and write, the stored values.
        """
        base = self.base
        suffix = BaseExperimentMasterSatellite.TWO_NODES_SUFFIX
        scanned = {self.scan_variable1, self.scan_variable2} - {None}
        before_scan = getattr(self, "_health_check_shared_values", {})
        forwarded = {}
        if self.EXPERIMENT_MODE == "two_nodes":
            for bare, shared in sorted(
                    base.two_node_shared_variable_map.items()):
                if shared in scanned and shared in before_scan:
                    forwarded[bare] = before_scan[shared]
                else:
                    forwarded[bare] = getattr(self, shared)
        node_names = base.node_variable_names[node]
        for target, value in self.authoritative_overrides.items():
            if target in scanned:
                continue
            if target in node_names:
                forwarded[target] = value
            elif (target in base.global_variable_names
                  and not target.endswith(suffix)):
                forwarded[target] = value
        return {name: value for name, value in forwarded.items()
                if name != "n_measurements"
                and not _MICROWAVE_CALIBRATION.match(name)}

    def _validate_forwarded_overrides(self, node):
        """Run the forwarded dictionary through the optimizer's own merge.

        The same merge_per_node_override_dictionaries call the optimizer's
        prepare() makes for *node*, with a single_node resolver for that node,
        so a name it would refuse stops this run here instead of crashing
        every check chain later.
        """
        text = format_override_dictionary(
            self._health_check_forwarded_overrides(node))
        resolver = BaseExperimentMasterSatellite(
            experiment=None, experiment_mode="single_node", which_node=node)
        fields = types.SimpleNamespace(**{
            f"override_ExperimentVariables_{each}":
                text if each == node else "{}"
            for each in self.VALID_NODES})
        try:
            return merge_per_node_override_dictionaries(
                experiment=fields,
                active_nodes=(node,),
                resolve_target=resolver.resolve_experiment_variable_target,
                eval_globals=globals(),
            )
        except ValueError as error:
            raise ValueError(
                f"run_health_check: the {node} optimizer would refuse this "
                f"run's settings ({text}): {error}") from error

    def _health_check_expid(self, node, scan_flag, forwarded_text):
        """The optimizer's expid for one check. Every key must be one of its
        arguments -- ARTIQ silently ignores any other, and a test runs this
        expid through the optimizer. What is not named keeps its default."""
        arguments = {
            "which_node": node,
            "parent_rid": int(self.scheduler.rid),
            "run_health_check_and_optimize": True,
            "target_fidelity": float(self.target_fidelity),
            "show_shared_applets": False,
            "n_measurements": int(self.health_check_n_measurements),
            "n_measurements_optimization_scan":
                int(self.health_check_optimization_n_measurements),
            "enable_fitting": True,
            "enable_geometric_frequency_scan": True,
            "enable_node_specific_adaptive_scan": True,
            "adaptive_half_range_linewidths":
                float(self.health_check_adaptive_half_range_linewidths),
            "adaptive_points_per_linewidth":
                float(self.health_check_adaptive_points_per_linewidth),
        }
        for flag in HEALTH_CHECK_OPTIMIZER_SCAN_FLAGS:
            arguments[flag] = flag == scan_flag
        for each in self.VALID_NODES:
            arguments[f"override_ExperimentVariables_{each}"] = (
                forwarded_text if each == node else "{}")
        return {
            "log_level": 30,
            "file": HEALTH_CHECK_OPTIMIZER_FILE,
            "class_name": HEALTH_CHECK_OPTIMIZER_CLASS,
            "arguments": arguments,
            "repo_rev": "N/A",
        }

    def _health_check_due(self, iteration):
        """Never before the first point: nothing has drifted yet."""
        return (bool(self.run_health_check) and iteration > 0
                and iteration % self._health_check_every_n == 0)

    def _initialize_health_check_records(self):
        """Present, if empty, in every file of a run with the check on.

        Also captures the shared two-node settings as they stand before the
        first scan point (overrides applied), for a scanned one's forwarding.
        """
        self._health_check_shared_values = {
            shared: getattr(self, shared) for shared in
            self.base.two_node_shared_variable_map.values()
            if hasattr(self, shared)}
        self._health_check_events = {key: [] for key in (
            "gvs_health_check_scan_point", "gvs_health_check_node",
            "gvs_health_check_transition", "gvs_health_check_fidelity",
            "gvs_health_check_resubmissions", "gvs_health_check_crashes",
            "gvs_health_check_center_before_Hz",
            "gvs_health_check_center_after_Hz",
            "gvs_health_check_time_start", "gvs_health_check_time_end")}
        self._health_check_chain_lines = []
        self._health_check_log = []
        keys = [(node, label) for node in self._health_check_nodes
                for label, _, _ in self._health_check_scans]
        # resubmissions after a failed check, and after a crashed chain
        self._health_check_totals = {key: 0 for key in keys}
        self._health_check_resubmitted_points = {key: [] for key in keys}
        self._health_check_crash_totals = {key: 0 for key in keys}
        self._health_check_crash_points = {key: [] for key in keys}
        # gvs_health_check_transition indexes this
        self.set_dataset("gvs_health_check_transition_labels",
                         ",".join(HEALTH_CHECK_LABELS),
                         broadcast=False, persist=False)
        self._save_health_check_records()
        self._publish_health_check_status(
            f"no check yet; the first is before scan point "
            f"{self._health_check_every_n}")
        self._create_health_check_status_applet()

    def _save_health_check_records(self):
        """LOCAL to this run's h5, never broadcast, and never named _Node1 or
        _Node2 (the analysis loaders strip those suffixes). Fresh lists every
        time, so an underflow rollback leaves them be."""
        for key, values in self._health_check_events.items():
            self.set_dataset(key, list(values), broadcast=False, persist=False)
        self.set_dataset("gvs_health_check_chain_rids",
                         "\n".join(self._health_check_chain_lines),
                         broadcast=False, persist=False)
        self.set_dataset("gvs_health_check_log",
                         "\n".join(self._health_check_log),
                         broadcast=False, persist=False)

    def _health_check_log_line(self, message, warn=False):
        line = f"{time.strftime('%H:%M:%S')} {message}"
        self._health_check_log.append(line)
        print("[health check] " + line)
        (logging.warning if warn else logging.info)("health check: %s", message)
        self.set_dataset("gvs_health_check_log",
                         "\n".join(self._health_check_log),
                         broadcast=False, persist=False)

    def _publish_health_check_status(self, current):
        """health_check_status, the status applet's text: what is happening
        now, and the resubmissions so far per node and transition, after a
        failed check and after a crashed chain counted apart."""
        totals = self._health_check_totals
        crash_totals = self._health_check_crash_totals
        lines = [
            f"GVS RID {self.scheduler.rid}: MW health check every "
            f"{self._health_check_every_n} scan point(s), "
            f"{', '.join(self._health_check_nodes)}: freq "
            f"{', '.join(label for label, _, _ in self._health_check_scans)}"
            f", target fidelity {float(self.target_fidelity):.2f}",
            current,
            "",
            "Resubmissions this run (failed check / crashed chain):",
        ]

        def at_points(points):
            return (f" (at scan point {', '.join(map(str, points))})"
                    if points else "")

        for (node, label), total in totals.items():
            points = self._health_check_resubmitted_points[(node, label)]
            crash_points = self._health_check_crash_points[(node, label)]
            lines.append(
                f"  {node} freq{label}: {total} failed{at_points(points)}, "
                f"{crash_totals[(node, label)]} crashed"
                f"{at_points(crash_points)}")
        lines.append(f"  total: {sum(totals.values())} failed, "
                     f"{sum(crash_totals.values())} crashed")
        self.set_dataset("health_check_status", "\n".join(lines),
                         broadcast=True, persist=False)

    def _create_health_check_status_applet(self):
        """The status applet, in a group no retirement pass touches. Never
        fatal, like _create_node_applets."""
        try:
            from applets_master_satellite import (
                HEALTH_CHECK_STATUS_APPLET_SPECS,
                create_shared_applets_for,
            )
            create_shared_applets_for(
                self, self.base,
                shared_specs=HEALTH_CHECK_STATUS_APPLET_SPECS)
        except Exception as error:  # noqa: BLE001 - convenience only
            logging.warning("could not create the health check applet: %s",
                            error)

    def _release_core(self):
        """Close the core device connection before every pause: the firmware
        drops a session when another connection opens, and the optimizer
        opens one. The next kernel call reconnects."""
        self.core.comm.close()

    def _raise_if_termination_requested(self):
        """Terminate pressed before a submission: stop before submitting.
        artiq_run's DummyScheduler has no check_termination."""
        check = getattr(self.scheduler, "check_termination", None)
        if check is not None and check():
            raise TerminationRequested

    def _health_check_runs(self):
        """{rid: status name} of every run the scheduler still holds that this
        GVS started: the optimizer class, with parent_rid = this RID. That
        includes the scan a failed check submits, which carries parent_rid
        on."""
        own = int(self.scheduler.rid)
        runs = {}
        for rid, run in self.scheduler.get_status().items():
            expid = run.get("expid") or {}
            if expid.get("class_name") != HEALTH_CHECK_OPTIMIZER_CLASS:
                continue
            try:
                parent = int(float(expid.get("arguments", {})["parent_rid"]))
            except (KeyError, TypeError, ValueError):
                continue
            if parent == own:
                runs[int(rid)] = run.get("status")
        return runs

    def _wait_for_health_check_chain(self, chain_rids):
        """Wait, inside this run, until the chain begun by chain_rids[0] is done.

        Appends every chain RID it sees. pause() hands the core over only to a
        PREPARED run of higher priority, so while the next one is still
        preparing this polls instead. Any queued run with a priority between
        this one's and 99 runs during the pause too (accepted, 2026-10-29).
        """
        first = chain_rids[0]
        while True:
            runs = self._health_check_runs()
            for rid in sorted(runs):
                # older RIDs belong to earlier checks still being analysed
                if rid >= first and rid not in chain_rids:
                    chain_rids.append(rid)
            if all(status in _HEALTH_CHECK_FINISHED
                   for rid, status in runs.items() if rid >= first):
                return
            # Terminate BEFORE check_pause(), while this run still holds the
            # pipeline: check_pause() is True for a Terminate too, but pause()
            # raises only when the master resumes this run -- after it has run
            # every prepared chain run, which with the lasers unlocked may
            # never finish (load-until-atom). check_termination() answers
            # without yielding, so the scan ends at once; the chain runs stay
            # queued and run after it (nothing is deleted).
            self._raise_if_termination_requested()
            if self.scheduler.check_pause():
                self._release_core()
                self.scheduler.pause()
            else:
                time.sleep(self.health_check_poll_interval)

    def _read_health_check_result(self, node, scan_flag, chain_rids):
        """(outcome, fidelity) of the chain, outcome "passed", "failed" or
        "crashed", from the optimizer's broadcast
        health_check_result_for_parent_<this RID> (archive=False: not this
        run's data).

        Keyed by this GVS's RID, never the latest health_check_result_*: a run
        queued between this run's priority and 99 -- a user's own optimizer,
        another GVS's chain -- runs while this one is paused and would
        overwrite those before this run looks. Only this GVS's chains write
        this key, one after another, so the result is this chain's if it is
        no older than the chain (a scan that finished before this run could
        see its RID is still this chain's) and is for this node and
        transition.

        "crashed": the chain ended with no result of its own -- none at all
        or an earlier attempt's (an exception in prepare or run), or only the
        CHECK's failure (chain_rids[0]): a failed check always submits a
        scan, and that scan (fitting on) always publishes its fit check, so
        the scan crashed (a FitError) or was never submitted. The fidelity is
        then NaN, or the failed check's.
        """
        result = self.get_dataset(
            f"health_check_result_for_parent_{int(self.scheduler.rid)}", None,
            archive=False)
        try:
            result_rid, result_node, result_scan, fidelity, passed = result
            result_rid = int(result_rid)
        except (TypeError, ValueError):
            return "crashed", float("nan")
        if (result_rid < chain_rids[0] or result_node != node
                or result_scan != scan_flag):
            return "crashed", float("nan")
        if result_rid not in chain_rids:
            chain_rids.append(result_rid)
        if bool(passed):
            return "passed", float(fidelity)
        if result_rid == chain_rids[0]:
            return "crashed", float(fidelity)
        return "failed", float(fidelity)

    def _health_check_crash_backoff(self):
        """Wait health_check_crash_backoff_s before resubmitting a crashed
        chain, so a crash that repeats at once cannot spin the queue. Asks
        check_termination() about once a second: a Terminate still ends the
        scan cleanly. A plain sleep, not pause(): nothing of this run's is
        queued to yield to."""
        remaining = float(self.health_check_crash_backoff_s)
        while True:
            self._raise_if_termination_requested()
            if remaining <= 0:
                return
            step = min(1.0, remaining)
            time.sleep(step)
            remaining -= step

    def _health_check_center_value(self, node, name):
        """The value this run uses; the master DB's for a node it does not."""
        if self._health_check_node_is_active(node):
            return float(getattr(self, name))
        value = self.get_dataset(name, None, archive=False)
        return float("nan") if value is None else float(value)

    def _run_one_health_check(self, iteration, node, label, scan_flag, center):
        """One (node, transition): check, re-optimize until passed, re-read."""
        what = f"{node} freq{label}"
        center_name = f"{center}_{node}"
        before = self._health_check_center_value(node, center_name)
        forwarded_text = format_override_dictionary(
            self._health_check_forwarded_overrides(node))
        expid = self._health_check_expid(node, scan_flag, forwarded_text)
        self._health_check_log_line(
            f"scan point {iteration}: {what}, forwarded "
            f"override_ExperimentVariables_{node} = {forwarded_text}")

        start = time.time()
        chain_rids = []
        resubmissions = 0
        crashes = 0
        while True:
            self._raise_if_termination_requested()
            attempt_rids = [int(self.scheduler.submit(
                pipeline_name="main", expid=copy.deepcopy(expid),
                priority=HEALTH_CHECK_PRIORITY, due_date=None, flush=False))]
            self._publish_health_check_status(
                f"{what}: checking since "
                f"{time.strftime('%H:%M:%S', time.localtime(start))}, "
                f"RID {attempt_rids[0]}, resubmissions this check "
                f"{resubmissions} ({crashes} after a crash)")
            self._wait_for_health_check_chain(attempt_rids)
            outcome, fidelity = self._read_health_check_result(
                node, scan_flag, attempt_rids)
            chain_rids += attempt_rids
            if outcome == "passed":
                break
            ### the SAME experiment again, for as long as it takes
            resubmissions += 1
            key = (node, label)
            if outcome == "failed":
                ### a genuine failure: at once, at 99, so it runs next
                self._health_check_totals[key] += 1
                points = self._health_check_resubmitted_points[key]
                if iteration not in points:
                    points.append(iteration)
                self._health_check_log_line(
                    f"scan point {iteration}: {what} did not pass (fidelity "
                    f"{fidelity:.3f}, RIDs {attempt_rids}); resubmitting now, "
                    f"resubmission {resubmissions} of this check", warn=True)
                continue
            ### no result of its own: the optimizer crashed (in prepare or
            ### run, a FitError, ...); wait first
            crashes += 1
            self._health_check_crash_totals[key] += 1
            points = self._health_check_crash_points[key]
            if iteration not in points:
                points.append(iteration)
            backoff = float(self.health_check_crash_backoff_s)
            self._health_check_log_line(
                f"scan point {iteration}: {what} chain RIDs {attempt_rids} "
                f"ended with no result of its own (crashed); resubmitting in "
                f"{backoff:.0f} s, resubmission {resubmissions} of this check "
                f"(crash {crashes})", warn=True)
            resume = time.localtime(time.time() + backoff)
            self._publish_health_check_status(
                f"{what}: chain RIDs {attempt_rids} crashed; resubmitting at "
                f"{time.strftime('%H:%M:%S', resume)}, resubmissions this "
                f"check {resubmissions} ({crashes} after a crash)")
            self._health_check_crash_backoff()

        ### only this center, only from the master DB: reload_experiment_
        ### variables would re-archive every variable, and
        ### _apply_run_wide_overrides re-write broadcast keys locally
        ### (deleting them from the master DB)
        after = self.get_dataset(center_name, None, archive=False)
        after = before if after is None else float(after)
        if self._health_check_node_is_active(node):
            setattr(self, center_name, after)
            self._refresh_runtime_variable_state()
        end = time.time()

        events = self._health_check_events
        events["gvs_health_check_scan_point"].append(iteration)
        events["gvs_health_check_node"].append(int(node[-1]))
        events["gvs_health_check_transition"].append(
            HEALTH_CHECK_LABELS.index(label))
        events["gvs_health_check_fidelity"].append(fidelity)
        events["gvs_health_check_resubmissions"].append(resubmissions)
        events["gvs_health_check_crashes"].append(crashes)
        events["gvs_health_check_center_before_Hz"].append(before)
        events["gvs_health_check_center_after_Hz"].append(after)
        events["gvs_health_check_time_start"].append(start)
        events["gvs_health_check_time_end"].append(end)
        self._health_check_chain_lines.append(
            f"scan point {iteration} {what}: "
            + ",".join(map(str, chain_rids)))
        self._save_health_check_records()
        self._health_check_log_line(
            f"scan point {iteration}: {what} passed, fidelity "
            f"{fidelity:.3f}, center {before:.0f} -> {after:.0f} Hz, "
            f"{resubmissions} resubmission(s) ({crashes} after a crash), "
            f"RIDs {chain_rids}")
        self._publish_health_check_status(
            f"{what}: passed at {time.strftime('%H:%M:%S')} (fidelity "
            f"{fidelity:.3f}) after {resubmissions} resubmission(s) "
            f"({crashes} after a crash)")

    def _run_health_check_event(self, iteration):
        """Every selected check before scan point *iteration*.

        Returns False when the user terminated the scan during it; the scan
        then ends (see _end_scan_on_health_check_termination).
        """
        if int(self.scheduler.rid) == 0:
            # artiq_run's DummyScheduler, and the offline compile sweep: no
            # master runs the optimizer, so there is nothing to wait for.
            self._health_check_log_line(
                f"scan point {iteration}: no ARTIQ master (RID 0), health "
                f"check skipped", warn=True)
            return True
        try:
            for node in self._health_check_nodes:
                for label, scan_flag, center in self._health_check_scans:
                    self._run_one_health_check(
                        iteration, node, label, scan_flag, center)
        except TerminationRequested:
            self._end_scan_on_health_check_termination(iteration)
            return False

        # single_node: the optimizer on the OTHER node left that node's
        # beams and coils as it finished; drive it safe again, as run() does
        # once at the start
        if self.EXPERIMENT_MODE == "single_node" and any(
                node != self.selected_node
                for node in self._health_check_nodes):
            self.force_other_node_off()
        self._restore_dashboard_after_health_check()
        self._publish_health_check_status(
            f"idle; the checks before scan point {iteration} finished at "
            f"{time.strftime('%H:%M:%S')}")
        self._save_scan_point_results()
        return True

    def _restore_dashboard_after_health_check(self):
        """The optimizer runs wrote their results, scan labels and
        n_measurements over this run's in the master DB, and their applets may
        have retired this run's. This run's h5 is unaffected."""
        republish_run_datasets(self._HasEnvironment__dataset_mgr)
        self._create_node_applets()
        self._create_health_check_status_applet()

    def _end_scan_on_health_check_termination(self, iteration):
        ### TERMINATE PRESSED DURING (OR JUST BEFORE) A HEALTH CHECK
        ### (confirmed by the user, 2026-10-29).
        ###
        ### TerminationRequested came from pause(), or from
        ### check_termination() while waiting, before a submission or in a
        ### crash backoff. DELETE NOTHING: every optimizer run of the chain
        ### stays in the queue and a running one finishes (the user: with
        ### the lasers unlocked nothing should run until they are relocked
        ### anyway, so leave them; "don't delete the existing queue"). Log
        ### which are left, save, and end the scan cleanly: run() returns
        ### without raising, so the h5 keeps every whole point.
        queued, running = [], []
        for rid, status in sorted(self._health_check_runs().items()):
            if status in _HEALTH_CHECK_NOT_STARTED:
                queued.append(rid)
            elif status not in _HEALTH_CHECK_FINISHED:
                running.append(rid)
        self._health_check_log_line(
            f"scan point {iteration}: terminated by the user during the "
            f"health check; nothing deleted: optimizer run(s) still queued "
            f"{queued}, running {running}; the scan ends before scan point "
            f"{iteration}", warn=True)
        self.set_dataset("gvs_health_check_terminated_at_scan_point",
                         iteration, broadcast=False, persist=False)
        self._publish_health_check_status(
            f"terminated by the user before scan point {iteration} at "
            f"{time.strftime('%H:%M:%S')}; optimizer RIDs left in the queue: "
            f"{queued + running}")
        republish_run_datasets(self._HasEnvironment__dataset_mgr)
        self._save_scan_point_results()

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
        if self.run_health_check:
            self._initialize_health_check_records()

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
            # before the point and outside its underflow retry, so a retried
            # point is never checked twice
            if self._health_check_due(iteration):
                if not self._run_health_check_event(iteration):
                    print("**************** General Variable Scan "
                          "master-satellite ENDED (terminated during the "
                          "health check) ****************")
                    return
            self._run_scan_point(variable1_value, variable2_value, iteration)
            iteration += 1

        print(
            "**************** General Variable Scan master-satellite DONE "
            "****************"
        )
