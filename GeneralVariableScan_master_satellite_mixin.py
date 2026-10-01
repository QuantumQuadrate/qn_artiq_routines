"""Shared scan implementation (mixins) for the master-satellite GVS files.

This module holds no public EnvExperiment class and is not an ARTIQ Explorer
experiment; the runnable experiments live in the mode-named GVS files.
"""

import inspect
import logging
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


class _GeneralVariableScanMasterSatelliteMixin(_DatasetRedirectMixin):
    """Shared implementation for the public master-satellite GVS files.

    This mixin deliberately does not inherit EnvExperiment, so ARTIQ Explorer
    cannot expose this common implementation as a third dashboard experiment.
    """

    VALID_NODES = ("Node1", "Node2")
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
        # One override dictionary per node, and deliberately no shared one:
        # both may stay populated at once, so switching selected_node is the
        # only edit needed to move a scan between nodes. Globals such as
        # n_measurements go in the running node's dictionary.
        for node in self.VALID_NODES:
            self.setattr_argument(
                f"override_ExperimentVariables_{node}",
                StringValue("{}"),
                tooltip=f"Overrides applied only when {node} runs. Otherwise "
                        f"ignored entirely -- not even evaluated, so it may "
                        f"name {node} variables that do not exist on the "
                        f"other node. Globals belong here too.",
            )
        # Optional per-point RTIOUnderflow retry. The three retry values
        # below take effect only when enable_Catch_UnderFlow is True.
        self.setattr_argument(
            "enable_Catch_UnderFlow",
            BooleanValue(False),
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
        self.setattr_argument(
            "skip_only_that_iteration_if_exhausted",
            BooleanValue(True),
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
        """Merge this run's per-node override dictionaries.

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
        if self.EXPERIMENT_MODE == "single_node":
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
        for target, value in self.authoritative_overrides.items():
            setattr(self, target, value)
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
        if self.EXPERIMENT_MODE == "single_node":
            # initialize_result_datasets() only covers the magnetometer
            # results. The reused atom-physics functions also need the full
            # single-node result surface (SPCM datasets, per-measurement
            # buffers and the host scalars they read in kernels), which the
            # microwave optimizer already sets up this way.
            self.base.initialize_single_node_result_state()

            # Scan labels, exactly as the standalone GeneralVariableScan
            # publishes them; the retention applet reads these for its axis.
            # Single-node only: Base publishes scan_var_dataset and friends
            # as part of the legacy-compatibility namespace, which two-node
            # mode does not have.
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

        Worth knowing about feedback RF in two-node mode: prepare_laser_
        stabilizer is deliberately skipped there, so such a run writes no
        p_AOM_A*_history_* at all and the plot shows the last single-node
        run's data rather than live values.

        Never fatal: the applets are a convenience, and the CCB is a
        dashboard-side service that artiq_run and the offline compile check
        do not provide meaningfully.
        """
        if not self.create_applets:
            return
        try:
            if self.EXPERIMENT_MODE == "single_node":
                from applets_master_satellite import create_applets_for

                created = create_applets_for(self, self.base)
            else:
                from applets_master_satellite import create_shared_applets_for

                created = create_shared_applets_for(self, self.base)
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

    def _report_underflow(self, message):
        logging.warning(message)
        print(message)

    def _underflow_backoff(self):
        """Pause on the host; the retried point resets the core via Base."""
        time.sleep(max(0.0, float(self.underflow_backoff_ms)) / 1000.0)

    def _run_scan_point(self, variable1_value, variable2_value, iteration):
        """Run one scan point, retrying only on RTIOUnderflow when enabled.

        With enable_Catch_UnderFlow off (the default) the point executes
        directly and any failure propagates unchanged. When it is on, only
        RTIOUnderflow is caught; a retried point repeats the full scan-point
        flow, including hardware initialization and its core reset. DRTIO,
        SPI, resolution, and all other failures always propagate.
        """
        if not self.enable_Catch_UnderFlow:
            self._execute_scan_point(
                variable1_value, variable2_value, iteration
            )
            return

        retries = 0
        while True:
            try:
                self._execute_scan_point(
                    variable1_value, variable2_value, iteration
                )
                return
            except RTIOUnderflow as error:
                retries += 1
                maximum_retries = int(self.underflow_max_retries)
                self._report_underflow(
                    f"RTIO underflow at iteration {iteration}, "
                    f"retry {retries}/{maximum_retries}: {error}"
                )
                if retries >= maximum_retries:
                    message = (
                        f"RTIO underflow at iteration {iteration} exceeded "
                        f"max retries ({maximum_retries})."
                    )
                    self._report_underflow(message)
                    if not self.skip_only_that_iteration_if_exhausted:
                        raise
                    return
                self._underflow_backoff()

    def run(self):
        self._initialize_run_state()

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
        for variable1_value in self.scan_sequence1:
            for variable2_value in self.scan_sequence2:
                self._run_scan_point(
                    variable1_value, variable2_value, iteration
                )
                iteration += 1

        print(
            "**************** General Variable Scan master-satellite DONE "
            "****************"
        )
