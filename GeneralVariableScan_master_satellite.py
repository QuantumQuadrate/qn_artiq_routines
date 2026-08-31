"""General variable scanning for the master-satellite execution stack."""

import inspect
import logging

import numpy as np
from numpy import array

from artiq.experiment import *

from utilities.BaseExperiment_master_satellite import (
    BaseExperimentMasterSatellite,
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
    import subroutines.experiment_functions_master_satellite as functions

    return build_experiment_function_registry(functions)


class GeneralVariableScanMasterSatellite(EnvExperiment):
    """GeneralVariableScan_master_satellite

    Scan authoritative variables on the master-satellite hardware stack.
    """

    VALID_MODES = ("single_node", "two_nodes")
    VALID_NODES = ("Node1", "Node2")
    LEGACY_NODE_NAMES = {"Node1": "alice", "Node2": "bob"}

    def build(self):
        self.setattr_argument(
            "experiment_mode", EnumerationValue(self.VALID_MODES)
        )
        self.setattr_argument(
            "selected_node", EnumerationValue(self.VALID_NODES)
        )

        selected_node = (
            self.selected_node
            if self.experiment_mode == "single_node"
            else None
        )
        self.base = BaseExperimentMasterSatellite(
            experiment=self,
            experiment_mode=self.experiment_mode,
            which_node=selected_node,
        )
        self.base.build()

        persistent_n_measurements = self.n_measurements
        self.setattr_argument(
            "n_measurements",
            NumberValue(
                persistent_n_measurements,
                ndecimals=0,
                step=1,
                type="int",
            ),
        )
        self.setattr_argument(
            "scan_variable1_name",
            StringValue(
                "f_FORT" if self.experiment_mode == "single_node"
                else "f_FORT_Node1"
            ),
        )
        self.setattr_argument(
            "scan_sequence1",
            StringValue(
                "np.array([self.f_FORT])"
                if self.experiment_mode == "single_node"
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
                f"{self.experiment_mode!r} mode."
            )
        self.setattr_argument(
            "experiment_function", EnumerationValue(function_names)
        )

        self._publish_legacy_node_compatibility()

    def _active_function_registry(self):
        if self.experiment_mode == "single_node":
            return build_single_node_function_registry()
        if self.experiment_mode == "two_nodes":
            return build_two_node_function_registry()
        raise ValueError(
            f"Unsupported experiment_mode {self.experiment_mode!r}; "
            "expected 'single_node' or 'two_nodes'."
        )

    def _publish_legacy_node_compatibility(self):
        if self.experiment_mode != "single_node":
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
        self.execution_n_measurements = self.n_measurements
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

    @rpc(flags={"async"})
    def append_to_dataset(self, name, value):
        """Redirect selected-node legacy magnetometer result names.

        Existing single-node magnetometer helpers use hard-coded unsuffixed
        strings. Only those known result names are translated; all other
        dataset appends retain normal EnvExperiment behavior.
        """
        resolved_name = self.base.resolve_result_dataset_name(name)
        super().append_to_dataset(resolved_name, value)

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

    def _run_scan_point(self, variable1_value, variable2_value, iteration):
        """Execute one scan point without rebuilding or preparing devices."""
        self._apply_scan_point(variable1_value, variable2_value)
        self.initialize_hardware()
        self.base.reset_result_state_for_scan_point()
        self.set_dataset("iteration", iteration, broadcast=True)
        self._selected_experiment_function(self)

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
