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
        self.setattr_argument(
            "override_ExperimentVariables", StringValue("{}")
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

        overrides = self._evaluate_expression(
            self.override_ExperimentVariables,
            "override_ExperimentVariables",
            {"self": self},
        )
        if not isinstance(overrides, dict):
            raise ValueError(
                "override_ExperimentVariables must evaluate to a dictionary."
            )
        self.authoritative_overrides = {
            self.base.resolve_experiment_variable_target(str(name)): value
            for name, value in overrides.items()
        }

        self.experiment_name = str(self.experiment_function)
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

    def _initialize_run_state(self):
        """Refresh queued values and initialize transient run results once."""
        if self.needs_experiment_variable_reload:
            self.base.reload_experiment_variables()
            # n_measurements is a run-local GUI value, not a persistent write.
            self.n_measurements = self.execution_n_measurements

        self._apply_run_wide_overrides()
        self.base.initialize_result_datasets()

    def _execute_scan_point(self, variable1_value, variable2_value, iteration):
        """Execute one scan point without rebuilding or preparing devices."""
        self._apply_scan_point(variable1_value, variable2_value)
        self.initialize_hardware()
        self.base.reset_result_state_for_scan_point()
        self.set_dataset("iteration", iteration, broadcast=True)
        self._selected_experiment_function(self)

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
