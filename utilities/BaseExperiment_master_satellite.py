"""Base hardware orchestration for the master-satellite architecture."""

import json
from pathlib import Path

import numpy as np
from artiq.experiment import kernel, delay, ms, rpc, us

from ExperimentVariables_master_satellite_Node1 import NODE1_VARIABLES
from ExperimentVariables_master_satellite_Node2 import NODE2_VARIABLES
from ExperimentVariables_master_satellite_global import (
    MASTER_SATELLITE_VARIABLES,
)
from utilities.DeviceAliases_master_satellite import (
    DeviceAliasesMasterSatellite,
)


class BaseExperimentMasterSatellite:
    """Expose and initialize master-satellite hardware for one or two nodes.

    Device resolution is delegated to :class:`DeviceAliasesMasterSatellite`.
    This class owns logical presentation, the fixed shared-SPCM policy, and the
    lifecycle of the one canonical core.
    """

    VALID_MODES = ("single_node", "two_nodes")
    VALID_NODES = ("Node1", "Node2")
    NODE_LEGACY_NAMES = {"Node1": "alice", "Node2": "bob"}
    SATELLITE_DESTINATION = 1

    DDS_ALIASES = (
        "dds_FORT",
        "dds_cooling_DP",
        "dds_D1_pumping_DP",
        "dds_MW_RF",
        "dds_microwaves",
        "GRIN1and2_dds",
        "dds_AOM_A1",
        "dds_AOM_A2",
        "dds_AOM_A3",
        "dds_AOM_A4",
        "dds_AOM_A5",
        "dds_AOM_A6",
    )

    CPLDS = ("urukul0_cpld", "urukul1_cpld", "urukul2_cpld")
    SAMPLERS = ("sampler0", "sampler1", "sampler2")
    ZOTINOS = ("zotino0",)
    TTLS = tuple(f"ttl{i}" for i in range(16))
    EDGE_COUNTERS = (
        "ttl0_counter",
        "ttl1_counter",
        "ttl2_counter",
        "ttl3_counter",
        "ttl8_counter",
        "ttl9_counter",
        "ttl10_counter",
        "ttl11_counter",
    )

    NODE_TTL_ALIASES = {
        "Node1": {
            "ttl_microwave_switch": "ttl4",
            "ttl_repump_switch": "ttl5",
            "ttl_exc0_switch": "ttl6",
            "ttl_pumping_repump_switch": "ttl7",
            "ttl_node2_input1": "ttl11",
            "ttl_node1_output1": "ttl12",
            "ttl_GRIN2_switch": "ttl13",
            "ttl_GRIN1_switch": "ttl14",
            "ttl_node1_output2": "ttl15",
        },
        "Node2": {
            "ttl_microwave_switch": "ttl4",
            "ttl_repump_switch": "ttl5",
            "ttl_exc0_switch": "ttl6",
            "ttl_pumping_repump_switch": "ttl7",
            "ttl_node1_input1": "ttl8",
            "ttl_node1_input2": "ttl9",
            "ttl_D1_pumping": "ttl12",
            "ttl_GRIN2_switch": "ttl13",
            "ttl_GRIN1_switch": "ttl14",
            "ttl_node2_output1": "ttl15",
        },
    }

    # Dummy aliases are retained only for legacy single-node presentation.
    COMPATIBILITY_DUMMY_TTL_ALIASES = {
        "Node1": {
            "ttl_node2_output1": "ttl15",
            "ttl_node1_input1": "ttl3",
            "ttl_node1_input2": "ttl3",
            "ttl_D1_pumping": "ttl3",
        },
        "Node2": {
            "ttl_node2_input1": "ttl11",
            "ttl_node1_output1": "ttl5",
            "ttl_node1_output2": "ttl5",
        },
    }

    SPCM_DEVICES = {
        "SPCM_H1": ("ttl0", "ttl0_counter"),
        "SPCM_V1": ("ttl1", "ttl1_counter"),
        "SPCM_H2": ("ttl8", "ttl8_counter"),
        "SPCM_V2": ("ttl9", "ttl9_counter"),
    }

    LEGACY_SPCM_ALIASES = {
        "ttl_SPCM0": "SPCM_H1",
        "ttl_SPCM1": "SPCM_V1",
        "ttl_SPCM0_OtherNode": "SPCM_H2",
        "ttl_SPCM1_OtherNode": "SPCM_V2",
    }

    SINGLE_NODE_WIRING_METADATA = {
        "Node1": {
            "coil_names": ["AZ bottom", "AZ top", "AX", "AY"],
            "AZ_bottom_Zotino_channel": 0,
            "AZ_top_Zotino_channel": 1,
            "AX_Zotino_channel": 13,
            "AY_Zotino_channel": 14,
            "coil_channels": [0, 1, 13, 14],
            "UV_trig_channel": [8],
            "Osc_trig_channel": [10],
            "FORT_MM_sampler_ch": 7,
            "GRIN1_sampler_ch": 4,
            "Magnetometer_X_ch": 1,
            "Magnetometer_Y_ch": 2,
            "Magnetometer_Z_ch": 3,
        },
        "Node2": {
            "coil_names": ["AZ bottom", "AZ top", "AX", "AY"],
            "AZ_bottom_Zotino_channel": 0,
            "AZ_top_Zotino_channel": 1,
            "AX_Zotino_channel": 2,
            "AY_Zotino_channel": 3,
            "coil_channels": [0, 1, 2, 3],
            "UV_trig_channel": [8],
            "Osc_trig_channel": [10],
            "FORT_MM_sampler_ch": 7,
            "GRIN1_sampler_ch": 4,
            "Magnetometer_X_ch": 1,
            "Magnetometer_Y_ch": 2,
            "Magnetometer_Z_ch": 3,
        },
    }

    MAGNETOMETER_RESULT_DATASETS = (
        "Magnetometer_Zero_X",
        "Magnetometer_Zero_Y",
        "Magnetometer_Zero_Z",
        "Magnetometer_OP_X",
        "Magnetometer_OP_Y",
        "Magnetometer_OP_Z",
        "Magnetometer_MOT_X",
        "Magnetometer_MOT_Y",
        "Magnetometer_MOT_Z",
    )

    # Legacy FORT-polarization helpers and routines address these transient
    # monitors and persistent waveplate calibrations by their standalone
    # names; single-node mode stores them node-suffixed.
    POLARIZATION_RESULT_DATASETS = (
        "FORT_MM_monitor",
        "FORT_APD_monitor",
        "HWP_angle",
        "QWP_angle",
        "best_852HWP_to_max",
        "best_852QWP_to_max",
        "best_852_power",
        "best_852_power_ref",
    )

    def __init__(
        self,
        experiment,
        experiment_mode=None,
        which_node=None,
        satellite_wait_attempts=100,
        satellite_ready_poll_interval=100 * ms,
    ):
        self.experiment = experiment
        self.experiment_mode = None
        self.which_node = None
        self.active_nodes = self.VALID_NODES
        self.node2_active = True
        self._execution_configured = False
        if experiment_mode is not None:
            self._set_execution_configuration(experiment_mode, which_node)
        elif which_node is not None:
            raise ValueError(
                "which_node cannot be supplied before experiment_mode is "
                "configured."
            )

        if satellite_wait_attempts < 1:
            raise ValueError("satellite_wait_attempts must be positive.")

        self.satellite_wait_attempts = int(satellite_wait_attempts)
        self.satellite_ready_poll_interval = satellite_ready_poll_interval

        self.node_resolvers = {}
        self.shared_spcm_resolver = None
        self.compatibility_variable_map = {}
        self.feedback_dataset_map = {}
        self.global_variable_names = frozenset(
            variable.name for variable in MASTER_SATELLITE_VARIABLES
        )
        self.node_variable_names = {
            node: frozenset(
                variable.name
                for variable in self._node_variable_definitions(node)
            )
            for node in self.VALID_NODES
        }
        self.unsuffixed_node_variable_names = frozenset(
            variable.name[:-len(f"_{node}")]
            for node in self.VALID_NODES
            for variable in self._node_variable_definitions(node)
        )
        self._built = False
        self._prepared = False

        self._cplds_node1 = []
        self._cplds_node2 = []
        self._samplers_node1 = []
        self._samplers_node2 = []
        self._zotinos_node1 = []
        self._zotinos_node2 = []
        self._dds_node1 = []
        self._dds_node2 = []
        self._dds_frequencies_node1 = []
        self._dds_frequencies_node2 = []
        self._dds_powers_node1 = []
        self._dds_powers_node2 = []
        self._ttl_inputs_node1 = []
        self._ttl_inputs_node2 = []
        self._ttl_outputs_node1 = []
        self._ttl_outputs_node2 = []
        self._ttl_safe_on_node1 = []
        self._ttl_safe_on_node2 = []
        self._ttl_safe_off_node1 = []
        self._ttl_safe_off_node2 = []
        self._spcm_inputs = []

    def _set_execution_configuration(self, experiment_mode, which_node):
        if experiment_mode not in self.VALID_MODES:
            raise ValueError(
                f"Unsupported master-satellite experiment_mode "
                f"{experiment_mode!r}; expected 'single_node' or 'two_nodes'."
            )
        if experiment_mode == "single_node":
            if which_node not in self.VALID_NODES:
                raise ValueError(
                    f"single_node mode requires which_node 'Node1' or "
                    f"'Node2', got {which_node!r}."
                )
            active_nodes = (which_node,)
        else:
            if which_node is not None:
                raise ValueError(
                    "which_node must be omitted in two_nodes mode."
                )
            active_nodes = self.VALID_NODES

        self.experiment_mode = experiment_mode
        self.which_node = which_node
        self.active_nodes = active_nodes
        self.node2_active = "Node2" in active_nodes
        self._execution_configured = True

    def configure_execution(self, experiment_mode, which_node=None):
        """Apply real submitted mode arguments after repository examination.

        An unconfigured Base binds a suffixed Node1/Node2 device superset in
        ``build()``.  This method then loads only the execution-relevant
        authoritative variables and publishes mode-specific presentation.
        """
        if self._prepared:
            raise RuntimeError("Cannot configure execution after prepare().")
        if self._execution_configured:
            if (
                self.experiment_mode != experiment_mode
                or self.which_node != which_node
            ):
                raise RuntimeError(
                    "Base execution mode is already configured as "
                    f"{self.experiment_mode!r}/{self.which_node!r}."
                )
            return

        self._set_execution_configuration(experiment_mode, which_node)
        if not self._built:
            return

        self._load_experiment_variables()
        if self.experiment_mode == "single_node":
            other_node = (
                "Node2" if self.which_node == "Node1" else "Node1"
            )
            self._deactivate_node_hardware_groups(other_node)
            self._publish_single_node_physical_presentations()
            self._bind_node_ttl_aliases(self.which_node)
            self._publish_shared_spcm_compatibility()
        self._publish_resolver_attributes()
        self._install_single_node_wiring_metadata()

    def _deactivate_node_hardware_groups(self, node):
        """Keep bound devices registered but exclude a node from lifecycle."""
        suffix = node.lower()
        for storage_name in (
            "_cplds",
            "_samplers",
            "_zotinos",
            "_ttl_inputs",
            "_ttl_outputs",
            "_ttl_safe_on",
            "_ttl_safe_off",
        ):
            setattr(self, f"{storage_name}_{suffix}", [])

    def _presentation_name(self, base_name, node):
        if self.experiment_mode == "single_node":
            return base_name
        return f"{base_name}_{node}"

    @staticmethod
    def _default_attribute_name(node, standalone_name):
        return f"{standalone_name}_{node}"

    def _make_resolver(self, node):
        return DeviceAliasesMasterSatellite(
            self.experiment,
            node,
            default_attribute_name_resolver=(
                lambda name, node=node: self._default_attribute_name(node, name)
            ),
        )

    def build(self):
        """Register and expose master-satellite physical devices."""
        if self._built:
            raise RuntimeError("BaseExperimentMasterSatellite.build() called twice.")

        if self._execution_configured:
            self._load_experiment_variables()

        for device_name in ("core", "core_dma", "scheduler"):
            self.experiment.setattr_device(device_name)

        nodes_to_bind = (
            self.active_nodes if self._execution_configured else self.VALID_NODES
        )
        for node in nodes_to_bind:
            resolver = self._make_resolver(node)
            self.node_resolvers[node] = resolver
            self._bind_node_physical_devices(node, resolver)
            self._bind_node_ttl_aliases(node)

        if "Node1" in self.node_resolvers:
            self.shared_spcm_resolver = self.node_resolvers["Node1"]
        else:
            # A Node2-only experiment still has access to shared master-local
            # SPCM hardware without constructing or initializing Node1 DDSes.
            self.shared_spcm_resolver = self._make_resolver("Node1")

        self._bind_shared_spcms()
        self._publish_resolver_attributes()
        self._install_single_node_wiring_metadata()
        self._built = True

    def _publish_single_node_physical_presentations(self):
        """Project the selected suffixed device superset to legacy names."""
        node = self.which_node
        for standalone_name in (
            self.CPLDS
            + self.SAMPLERS
            + self.ZOTINOS
            + self.TTLS
            + self.EDGE_COUNTERS
        ):
            setattr(
                self.experiment,
                standalone_name,
                getattr(self.experiment, f"{standalone_name}_{node}"),
            )

    def _install_single_node_wiring_metadata(self):
        """Expose fixed legacy wiring metadata for the selected node."""
        if self.experiment_mode != "single_node":
            return

        for name, value in self.SINGLE_NODE_WIRING_METADATA[
            self.which_node
        ].items():
            # Copy lists so experiment code cannot mutate the class constant.
            setattr(
                self.experiment,
                name,
                list(value) if isinstance(value, list) else value,
            )

        self.experiment.measurements_progress = "measurements_progress"

    def initialize_result_datasets(self):
        """Create the minimal transient single-node magnetometer results."""
        if self.experiment_mode != "single_node":
            return
        self.reset_result_state_for_scan_point()

    def resolve_result_dataset_name(self, name):
        """Resolve hard-coded legacy result names for the selected node."""
        if self.experiment_mode == "single_node":
            if (
                name in self.MAGNETOMETER_RESULT_DATASETS
                or name in self.POLARIZATION_RESULT_DATASETS
            ):
                return f"{name}_{self.which_node}"
            resolved_feedback_name = self.feedback_dataset_map.get(name)
            if resolved_feedback_name is not None:
                return resolved_feedback_name
        return name

    def reset_result_state_for_scan_point(self):
        """Clear transient magnetometer results for one scan point.

        This method deliberately does not read ExperimentVariables, touch the
        core, bind devices, or reset any variable-dependent cached state.
        """
        if self.experiment_mode != "single_node":
            return

        self.experiment.set_dataset(
            self.experiment.measurements_progress,
            0.0,
            broadcast=True,
            persist=False,
        )
        for legacy_dataset_name in self.MAGNETOMETER_RESULT_DATASETS:
            dataset_name = self.resolve_result_dataset_name(
                legacy_dataset_name
            )
            self.experiment.set_dataset(
                dataset_name,
                [0.0],
                broadcast=True,
                persist=False,
            )

    @staticmethod
    def _node_variable_definitions(node):
        return {
            "Node1": NODE1_VARIABLES,
            "Node2": NODE2_VARIABLES,
        }[node]

    def _load_experiment_variables(self):
        """Load authoritative datasets without creating or modifying them."""
        missing_by_owner = {}

        for node in self.active_nodes:
            missing = self._load_variable_definitions(
                self._node_variable_definitions(node)
            )
            if missing:
                missing_by_owner[
                    f"ExperimentVariables_master_satellite_{node}.py"
                ] = missing

        missing_globals = self._load_variable_definitions(
            MASTER_SATELLITE_VARIABLES
        )
        if missing_globals:
            missing_by_owner[
                "ExperimentVariables_master_satellite_global.py"
            ] = missing_globals

        if missing_by_owner:
            self._raise_missing_datasets(missing_by_owner)

        if self.experiment_mode == "single_node":
            suffix = f"_{self.which_node}"
            self.compatibility_variable_map = {
                variable.name[:-len(suffix)]: variable.name
                for variable in self._node_variable_definitions(
                    self.which_node
                )
            }
            self.refresh_compatibility_variables()

    def _load_variable_definitions(self, definitions):
        missing = []
        for variable in definitions:
            try:
                value = self.experiment.get_dataset(variable.name)
            except KeyError:
                missing.append(variable.name)
                continue
            setattr(self.experiment, variable.name, value)
        return missing

    def refresh_compatibility_variables(self):
        """Refresh single-node legacy attributes from authoritative values.

        This method only copies in-memory experiment attributes. It never
        reads, creates, or writes persistent datasets.
        """
        if self.experiment_mode != "single_node":
            return

        for compatibility_name, authoritative_name in (
            self.compatibility_variable_map.items()
        ):
            try:
                value = getattr(self.experiment, authoritative_name)
            except AttributeError as error:
                raise RuntimeError(
                    "Missing authoritative master-satellite experiment "
                    f"attribute {authoritative_name!r}."
                ) from error
            setattr(self.experiment, compatibility_name, value)

    def resolve_experiment_variable_target(self, name):
        """Resolve a GVS-facing name to its authoritative attribute name."""
        if not isinstance(name, str) or not name:
            raise ValueError(
                "Experiment-variable target name must be a non-empty string."
            )

        if name in self.global_variable_names:
            return name

        if self.experiment_mode == "single_node":
            selected_names = self.node_variable_names[self.which_node]
            if name in selected_names:
                return name

            other_node = "Node2" if self.which_node == "Node1" else "Node1"
            if name in self.node_variable_names[other_node]:
                raise ValueError(
                    f"Experiment-variable target {name!r} belongs to "
                    f"{other_node}, but single_node mode selected "
                    f"{self.which_node}."
                )

            try:
                return self.compatibility_variable_map[name]
            except KeyError:
                raise ValueError(
                    f"Unknown experiment-variable target {name!r} for "
                    f"single_node {self.which_node}."
                ) from None

        for node in self.VALID_NODES:
            if name in self.node_variable_names[node]:
                return name

        if name in self.unsuffixed_node_variable_names:
            raise ValueError(
                f"Ambiguous node-specific experiment-variable target "
                f"{name!r} in two_nodes mode; use {name + '_Node1'!r} or "
                f"{name + '_Node2'!r}."
            )

        raise ValueError(
            f"Unknown experiment-variable target {name!r} in two_nodes mode."
        )

    def refresh_variable_dependent_state(self):
        """Refresh cached DDS defaults from current authoritative attributes."""
        if not self._prepared:
            raise RuntimeError(
                "prepare() must be called before refreshing variable-dependent "
                "state."
            )

        for node in self.active_nodes:
            resolver = self.node_resolvers[node]
            frequencies = []
            powers = []
            for binding in resolver.dds_bindings:
                frequency_name = binding["frequency_attribute"]
                power_name = binding["power_attribute"]
                try:
                    frequencies.append(
                        getattr(self.experiment, frequency_name)
                    )
                    powers.append(getattr(self.experiment, power_name))
                except AttributeError as error:
                    raise RuntimeError(
                        f"Missing authoritative DDS variable while refreshing "
                        f"{node} alias {binding['logical_alias']!r}: {error}"
                    ) from error

            resolver.dds_frequencies[:] = frequencies
            resolver.dds_powers[:] = powers

    def reload_experiment_variables(self):
        """Reload active authoritative/global datasets without side effects."""
        if not self._built or not self._prepared:
            raise RuntimeError(
                "build() and prepare() must complete before reloading "
                "master-satellite experiment variables."
            )

        values = {}
        missing_by_owner = {}
        for node in self.active_nodes:
            node_values, missing = self._read_variable_definitions(
                self._node_variable_definitions(node)
            )
            values.update(node_values)
            if missing:
                missing_by_owner[
                    f"ExperimentVariables_master_satellite_{node}.py"
                ] = missing

        global_values, missing_globals = self._read_variable_definitions(
            MASTER_SATELLITE_VARIABLES
        )
        values.update(global_values)
        if missing_globals:
            missing_by_owner[
                "ExperimentVariables_master_satellite_global.py"
            ] = missing_globals

        if missing_by_owner:
            self._raise_missing_datasets(missing_by_owner)

        for name, value in values.items():
            setattr(self.experiment, name, value)

        self.refresh_compatibility_variables()
        self.refresh_variable_dependent_state()

    def _read_variable_definitions(self, definitions):
        values = {}
        missing = []
        for variable in definitions:
            try:
                values[variable.name] = self.experiment.get_dataset(
                    variable.name
                )
            except KeyError:
                missing.append(variable.name)
        return values, missing

    @staticmethod
    def _raise_missing_datasets(missing_by_owner):
        details = []
        for owner, names in missing_by_owner.items():
            details.append(f"{owner}: {', '.join(names)}")
        raise RuntimeError(
            "Missing required master-satellite persistent datasets. "
            "Run the indicated ExperimentVariables initializer(s): "
            + "; ".join(details)
        )

    def _bind_node_physical_devices(self, node, resolver):
        groups = (
            (self.CPLDS, "_cplds"),
            (self.SAMPLERS, "_samplers"),
            (self.ZOTINOS, "_zotinos"),
            (self.TTLS, None),
            (self.EDGE_COUNTERS, None),
        )
        for standalone_names, storage_prefix in groups:
            bound_devices = []
            for standalone_name in standalone_names:
                presentation_name = self._presentation_name(
                    standalone_name, node
                )
                device = resolver.bind_physical_device(
                    standalone_name, presentation_name
                )
                bound_devices.append(device)
            if storage_prefix is not None:
                setattr(
                    self,
                    f"{storage_prefix}_{node.lower()}",
                    bound_devices,
                )

        input_names = {
            "Node1": ("ttl11",),
            "Node2": ("ttl0", "ttl1", "ttl2", "ttl3", "ttl8", "ttl9"),
        }[node]
        output_names = {
            "Node1": (
                "ttl4", "ttl5", "ttl6", "ttl7", "ttl12", "ttl13", "ttl14", "ttl15"
            ),
            "Node2": (
                "ttl4", "ttl5", "ttl6", "ttl7", "ttl12", "ttl13", "ttl14", "ttl15"
            ),
        }[node]
        safe_on_names = {
            "Node1": ("ttl4", "ttl5", "ttl6", "ttl7", "ttl13", "ttl14"),
            "Node2": ("ttl4", "ttl5", "ttl6", "ttl7", "ttl12", "ttl13", "ttl14"),
        }[node]
        safe_off_names = {
            "Node1": ("ttl12", "ttl15"),
            "Node2": ("ttl15",),
        }[node]

        def references(names):
            return [
                getattr(self.experiment, self._presentation_name(name, node))
                for name in names
            ]

        suffix = node.lower()
        setattr(self, f"_ttl_inputs_{suffix}", references(input_names))
        setattr(self, f"_ttl_outputs_{suffix}", references(output_names))
        setattr(self, f"_ttl_safe_on_{suffix}", references(safe_on_names))
        setattr(self, f"_ttl_safe_off_{suffix}", references(safe_off_names))

    def _bind_node_ttl_aliases(self, node):
        mappings = dict(self.NODE_TTL_ALIASES[node])
        if self.experiment_mode == "single_node":
            mappings.update(self.COMPATIBILITY_DUMMY_TTL_ALIASES[node])

        for logical_name, standalone_name in mappings.items():
            physical_attribute = self._presentation_name(standalone_name, node)
            logical_attribute = self._presentation_name(logical_name, node)
            setattr(
                self.experiment,
                logical_attribute,
                getattr(self.experiment, physical_attribute),
            )

    def _bind_shared_spcms(self):
        # Raw TTL names retain ordinary per-node translation. Canonical SPCMs
        # are a separate detector-identity policy and ALWAYS resolve through
        # the Node1/master physical map. In particular, Node2 ttl0/1/8/9 are
        # not active SPCMs, and deprecated code using those raw names is not
        # rerouted here.
        for canonical_name, (ttl_name, counter_name) in self.SPCM_DEVICES.items():
            ttl_device = self.shared_spcm_resolver.bind_physical_device(
                ttl_name, canonical_name
            )
            counter_device = self.shared_spcm_resolver.bind_physical_device(
                counter_name, f"{canonical_name}_counter"
            )
            self._spcm_inputs.append(ttl_device)

        if self.experiment_mode == "single_node":
            self._publish_shared_spcm_compatibility()

    def _publish_shared_spcm_compatibility(self):
        """Publish fixed legacy SPCM aliases after mode configuration."""
        for legacy_name, canonical_name in self.LEGACY_SPCM_ALIASES.items():
            setattr(
                self.experiment,
                legacy_name,
                getattr(self.experiment, canonical_name),
            )
            setattr(
                self.experiment,
                f"{legacy_name}_counter",
                getattr(self.experiment, f"{canonical_name}_counter"),
            )

        if self.which_node == "Node2":
            # setattr_device("ttl0"), etc. necessarily creates the low-level
            # master attribute while canonical SPCMs are captured. Restore the
            # ordinary Node2 raw-TTL presentation afterward; the canonical and
            # legacy SPCM aliases above retain their master-local references.
            node2_resolver = self.node_resolvers["Node2"]
            for ttl_name, counter_name in self.SPCM_DEVICES.values():
                for standalone_name in (ttl_name, counter_name):
                    unified_name = node2_resolver.resolve_physical_name(
                        standalone_name
                    )
                    setattr(
                        self.experiment,
                        standalone_name,
                        getattr(self.experiment, unified_name),
                    )

    def _publish_resolver_attributes(self):
        if self.experiment_mode == "single_node":
            resolver = self.node_resolvers[self.which_node]
            self.experiment.named_devices = resolver
        else:
            self.experiment.named_devices_Node1 = self.node_resolvers["Node1"]
            self.experiment.named_devices_Node2 = self.node_resolvers["Node2"]

    def prepare(self):
        """Bind DDS aliases after suffixed authoritative variables exist."""
        if not self._built:
            raise RuntimeError("build() must be called before prepare().")
        if not self._execution_configured:
            raise RuntimeError(
                "configure_execution() must be called with submitted mode "
                "arguments before prepare()."
            )
        if self._prepared:
            raise RuntimeError("BaseExperimentMasterSatellite.prepare() called twice.")

        for node in self.active_nodes:
            resolver = self.node_resolvers[node]
            self.refresh_compatibility_variables()
            for logical_alias in self.DDS_ALIASES:
                resolver.bind_dds(
                    logical_alias,
                    self._presentation_name(logical_alias, node),
                )

            suffix = node.lower()
            setattr(self, f"_dds_{suffix}", resolver.dds_list)
            setattr(
                self,
                f"_dds_frequencies_{suffix}",
                resolver.dds_frequencies,
            )
            setattr(self, f"_dds_powers_{suffix}", resolver.dds_powers)

            all_dds_attribute = self._presentation_name("all_dds_channels", node)
            setattr(self.experiment, all_dds_attribute, list(resolver.dds_list))

        self._prepared = True

    def _feedback_channel_dataset_names(self):
        """List the selected node's feedback dataset names from its config."""
        config_path = (
            Path(__file__).resolve().parent
            / "config"
            / self.NODE_LEGACY_NAMES[self.which_node]
            / "feedback_channels.json"
        )
        with config_path.open() as config_file:
            stabilizer_dict = json.load(config_file)

        dataset_names = []
        for feedback_channels in stabilizer_dict.values():
            for channel_parameters in feedback_channels.values():
                dataset_names.append(channel_parameters["dataset"])
                dataset_names.append(channel_parameters["power_dataset"])
                dataset_names.append(
                    channel_parameters["power_dataset"] + "_history"
                )
        return dataset_names

    def prepare_laser_stabilizer(self, stabilizer_factory=None):
        """Build the selected node's AOMPowerStabilizer with suffixed persistence.

        subroutines/aom_feedback.py is unchanged from the standalone stack:
        every dataset it touches goes through the experiment object, so this
        method installs a name map that routes those reads and writes to the
        node-suffixed authoritative datasets before construction. The
        stabilizer itself sees the ordinary projected legacy namespace
        (samplers, DDS aliases, set points, dds_defaults, which_node).
        """
        if not self._prepared:
            raise RuntimeError(
                "prepare() must complete before preparing the laser "
                "stabilizer."
            )
        if self.experiment_mode != "single_node":
            raise RuntimeError(
                "Master-satellite laser feedback is single-node only for now."
            )

        suffix = f"_{self.which_node}"
        self.feedback_dataset_map = {
            name: f"{name}{suffix}"
            for name in self._feedback_channel_dataset_names()
        }
        self.feedback_dataset_map["feedbackchannels"] = (
            f"feedbackchannels{suffix}"
        )

        if stabilizer_factory is None:
            from subroutines.aom_feedback import AOMPowerStabilizer

            stabilizer_factory = AOMPowerStabilizer

        experiment = self.experiment
        # Present the resolver's DDS default-variable names the same way the
        # standalone DeviceAliases does for the stabilizer.
        experiment.dds_defaults = (
            self.node_resolvers[self.which_node].dds_defaults
        )
        fast_feedback_dds_names = eval(experiment.fast_feedback_dds_list)
        experiment.laser_stabilizer = stabilizer_factory(
            experiment=experiment,
            dds_names=fast_feedback_dds_names,
            iterations=experiment.aom_feedback_iterations,
            averages=experiment.aom_feedback_averages,
            leave_AOMs_on=False,
            leave_MOT_AOMs_on=True,
        )

        channels = experiment.laser_stabilizer.all_channels
        experiment.set_dataset(
            "feedbackchannels",
            [
                self.resolve_result_dataset_name(channel.dB_dataset)
                for channel in channels
            ],
            broadcast=True,
            persist=True,
        )
        experiment.initial_RF_dB_values = np.zeros(len(channels))
        for channel_index, channel in enumerate(channels):
            experiment.initial_RF_dB_values[channel_index] = (
                experiment.get_dataset(channel.dB_dataset, archive=False)
            )
            try:
                experiment.append_to_dataset(
                    channel.dB_history_dataset,
                    float(experiment.initial_RF_dB_values[channel_index]),
                )
            except KeyError:
                experiment.set_dataset(
                    channel.dB_history_dataset,
                    [float(experiment.initial_RF_dB_values[channel_index])],
                    broadcast=True,
                )
        return experiment.laser_stabilizer

    @kernel
    def _wait_for_satellite(self):
        ready = False
        for _ in range(self.satellite_wait_attempts):
            if self.experiment.core.get_rtio_destination_status(
                self.SATELLITE_DESTINATION
            ):
                ready = True
                break
            delay(self.satellite_ready_poll_interval)
        if not ready:
            raise RuntimeError("DRTIO destination 1 did not become available.")
        self.experiment.core.break_realtime()

    @kernel
    def initialize_hardware(
        self,
        turn_off_dds_channels=True,
        turn_off_zotinos=True,
        initialize_spcms=True,
        reset_core=True,
    ):
        """Initialize active hardware with ordered access and optional reset."""
        if not self._prepared:
            raise RuntimeError("prepare() must be called before initialize_hardware().")

        # This is the only core reset owned by the master-satellite base.
        # Node-local manual controls may omit it to preserve outputs already
        # established on the other node through the shared RTIO fabric.
        if reset_core:
            self.experiment.core.reset()

        if self.node2_active:
            self._wait_for_satellite()
        else:
            self.experiment.core.break_realtime()

        if initialize_spcms:
            for spcm in self._spcm_inputs:
                spcm.input()
        self.experiment.core.break_realtime()

        self._initialize_cplds(self._cplds_node1)
        self._initialize_cplds(self._cplds_node2)

        self._initialize_dds_group(
            self._dds_node1,
            self._dds_frequencies_node1,
            self._dds_powers_node1,
            turn_off_dds_channels,
        )
        self._initialize_dds_group(
            self._dds_node2,
            self._dds_frequencies_node2,
            self._dds_powers_node2,
            turn_off_dds_channels,
        )

        self._configure_ttl_group(
            self._ttl_inputs_node1,
            self._ttl_outputs_node1,
            self._ttl_safe_on_node1,
            self._ttl_safe_off_node1,
        )
        self._configure_ttl_group(
            self._ttl_inputs_node2,
            self._ttl_outputs_node2,
            self._ttl_safe_on_node2,
            self._ttl_safe_off_node2,
        )

        self._initialize_sampler_group(self._samplers_node1, True)
        self._initialize_sampler_group(self._samplers_node2, False)

        self._initialize_zotino_group(
            self._zotinos_node1, turn_off_zotinos
        )
        self._initialize_zotino_group(
            self._zotinos_node2, turn_off_zotinos
        )
        self.experiment.core.break_realtime()

    @kernel
    def _initialize_cplds(self, cplds):
        for cpld in cplds:
            self.experiment.core.break_realtime()
            cpld.init()
            delay(1 * ms)

    @kernel
    def _initialize_dds_group(
        self, devices, frequencies, powers, turn_off_outputs
    ):
        for i in range(len(devices)):
            self.experiment.core.break_realtime()
            devices[i].init()
            if turn_off_outputs:
                # Establish the safe switch state before changing attenuation
                # or profile data.  This is important for deterministic manual
                # hardware bring-up after a previous experiment left an output
                # enabled.
                devices[i].sw.off()
            devices[i].set_att(0.0)
            amplitude = (2 * 50 * 10 ** (powers[i] / 10 - 3)) ** 0.5
            devices[i].set(
                frequency=frequencies[i], amplitude=amplitude
            )
            delay(1 * ms)

    @kernel
    def _configure_ttl_group(self, inputs, outputs, safe_on, safe_off):
        for ttl in inputs:
            ttl.input()
        for ttl in outputs:
            ttl.output()
        for ttl in safe_on:
            ttl.on()
            delay(1 * ms)
        for ttl in safe_off:
            ttl.off()
            delay(1 * ms)
        self.experiment.core.break_realtime()

    @kernel
    def _initialize_sampler_group(self, samplers, configure_node1_gains):
        for sampler in samplers:
            self.experiment.core.break_realtime()
            sampler.init()
        if configure_node1_gains and len(samplers) >= 2:
            for channel in range(8):
                samplers[0].set_gain_mu(channel, 8)
                delay(100 * us)
            for channel in range(8):
                samplers[1].set_gain_mu(channel, 0)
                delay(100 * us)

    @kernel
    def _initialize_zotino_group(self, zotinos, turn_off_outputs):
        for zotino in zotinos:
            self.experiment.core.break_realtime()
            zotino.init()
            if turn_off_outputs:
                for channel in range(16):
                    zotino.write_dac(channel, 0.0)
                    zotino.load()
                    delay(1 * ms)


class _DatasetRedirectMixin:
    """Route legacy hard-coded dataset names through Base resolution.

    Reused single-node code addresses magnetometer results and laser-feedback
    datasets by their standalone names. These overrides translate exactly the
    names Base knows about (see resolve_result_dataset_name) and pass every
    other name through unchanged. This mixin deliberately does not inherit
    EnvExperiment and expects the experiment to hold its Base as self.base.
    """

    def set_dataset(self, key, value, **kwargs):
        return super().set_dataset(
            self.base.resolve_result_dataset_name(key), value, **kwargs
        )

    def get_dataset(self, key, *args, **kwargs):
        return super().get_dataset(
            self.base.resolve_result_dataset_name(key), *args, **kwargs
        )

    @rpc(flags={"async"})
    def append_to_dataset(self, key, value):
        super().append_to_dataset(
            self.base.resolve_result_dataset_name(key), value
        )
