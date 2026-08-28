"""Device resolution for the master-satellite hardware architecture.

This module translates a node's standalone physical device names through the
explicit master-satellite mapping and binds the resulting device objects to
caller-selected experiment attributes.  It deliberately does not manage the
core, satellite readiness, or hardware initialization order.
"""

import json
from pathlib import Path


_CONFIG_DIRECTORY = Path(__file__).resolve().parent / "config"
_NODE_CONFIG_DIRECTORIES = {
    "Node1": "alice",
    "Node2": "bob",
}


class DeviceResolutionError(RuntimeError):
    """Raised when a required logical, physical, or experiment binding fails."""


class DeviceAliasesMasterSatellite:
    """Resolve and bind devices for one master-satellite node context.

    ``default_attribute_name_resolver`` converts the standalone variable names
    in ``DDS_DEFAULTS`` into the node-specific experiment attributes that are
    the source of truth.  For example, a Node2 resolver can convert ``f_FORT``
    to ``f_FORT_Node2`` even when the DDS itself is presented as the legacy
    compatibility attribute ``dds_FORT``.
    """

    def __init__(
        self,
        experiment,
        node,
        default_attribute_name_resolver=None,
        logical_config=None,
        physical_mapping=None,
    ):
        if node not in _NODE_CONFIG_DIRECTORIES:
            raise ValueError(
                f"Unsupported master-satellite node {node!r}; "
                "expected 'Node1' or 'Node2'."
            )

        self.experiment = experiment
        self.node = node
        self.default_attribute_name_resolver = default_attribute_name_resolver

        if logical_config is None:
            logical_config = self._load_json(
                _CONFIG_DIRECTORY
                / _NODE_CONFIG_DIRECTORIES[node]
                / "device_aliases.json"
            )
        if physical_mapping is None:
            physical_mapping = self._load_json(
                _CONFIG_DIRECTORY / "master_satellite" / "device_aliases.json"
            )

        try:
            self.alias_map = logical_config["ALIAS_MAP"]
            self.dds_defaults = logical_config["DDS_DEFAULTS"]
            self.physical_map = physical_mapping[node]
        except KeyError as error:
            raise DeviceResolutionError(
                f"Missing required configuration section {error.args[0]!r} "
                f"for {node}."
            ) from error

        # These collections belong to this resolver instance, so Node1 and
        # Node2 bookkeeping cannot overwrite one another.
        self.dds_list = []
        self.dds_powers = []
        self.dds_frequencies = []
        self.dds_names_aliases = []
        self.dds_bindings = []

    @staticmethod
    def _load_json(path):
        try:
            with path.open() as config_file:
                return json.load(config_file)
        except (OSError, json.JSONDecodeError) as error:
            raise DeviceResolutionError(
                f"Unable to load device configuration {str(path)!r}: {error}"
            ) from error

    def resolve_physical_name(self, standalone_name):
        """Return the exact unified device_db name for a standalone name."""
        try:
            return self.physical_map[standalone_name]
        except KeyError as error:
            raise DeviceResolutionError(
                f"No master-satellite physical mapping for {self.node} "
                f"standalone device {standalone_name!r}."
            ) from error

    def bind_physical_device(self, standalone_name, experiment_attribute_name):
        """Resolve a physical device and expose it under a caller-chosen name."""
        unified_name = self.resolve_physical_name(standalone_name)

        try:
            self.experiment.setattr_device(unified_name)
            device = getattr(self.experiment, unified_name)
        except Exception as error:
            raise DeviceResolutionError(
                f"Unable to bind {self.node} standalone device "
                f"{standalone_name!r} to unified device {unified_name!r}: {error}"
            ) from error

        setattr(self.experiment, experiment_attribute_name, device)
        return device

    def resolve_dds_physical_name(self, logical_alias):
        """Resolve a logical DDS alias to standalone and unified names."""
        try:
            standalone_name = self.alias_map[logical_alias]
        except KeyError as error:
            raise DeviceResolutionError(
                f"DDS alias {logical_alias!r} is not defined for {self.node}."
            ) from error

        unified_name = self.resolve_physical_name(standalone_name)
        return standalone_name, unified_name

    def resolve_dds_default_attributes(
        self, logical_alias, default_attribute_name_resolver=None
    ):
        """Resolve a DDS alias's frequency and power source attributes.

        The supplied resolver maps legacy names from ``DDS_DEFAULTS`` to the
        globally unique master-satellite attributes.  Raw unsuffixed names are
        never assumed to be persistent master-satellite sources of truth.
        """
        try:
            defaults = self.dds_defaults[logical_alias]
            frequency_name = defaults["frequency"]
            power_name = defaults["power"]
        except KeyError as error:
            raise DeviceResolutionError(
                f"Incomplete DDS defaults for alias {logical_alias!r} on "
                f"{self.node}; missing {error.args[0]!r}."
            ) from error

        name_resolver = (
            default_attribute_name_resolver
            if default_attribute_name_resolver is not None
            else self.default_attribute_name_resolver
        )
        if name_resolver is None:
            raise DeviceResolutionError(
                f"No default-variable attribute resolver supplied for "
                f"{self.node} DDS alias {logical_alias!r}."
            )

        return name_resolver(frequency_name), name_resolver(power_name)

    def bind_dds(
        self,
        logical_alias,
        experiment_attribute_name,
        default_attribute_name_resolver=None,
    ):
        """Resolve, bind, and record one logical DDS device."""
        standalone_name, unified_name = self.resolve_dds_physical_name(
            logical_alias
        )
        frequency_attribute, power_attribute = (
            self.resolve_dds_default_attributes(
                logical_alias, default_attribute_name_resolver
            )
        )

        try:
            frequency = getattr(self.experiment, frequency_attribute)
            power = getattr(self.experiment, power_attribute)
        except AttributeError as error:
            raise DeviceResolutionError(
                f"Missing node-specific DDS default attribute for {self.node} "
                f"alias {logical_alias!r}: {error}"
            ) from error

        device = self.bind_physical_device(
            standalone_name, experiment_attribute_name
        )

        self.dds_list.append(device)
        self.dds_frequencies.append(frequency)
        self.dds_powers.append(power)
        self.dds_names_aliases.append(
            (experiment_attribute_name, unified_name)
        )
        self.dds_bindings.append(
            {
                "logical_alias": logical_alias,
                "experiment_attribute": experiment_attribute_name,
                "standalone_name": standalone_name,
                "unified_name": unified_name,
                "frequency_attribute": frequency_attribute,
                "power_attribute": power_attribute,
            }
        )
        return device
