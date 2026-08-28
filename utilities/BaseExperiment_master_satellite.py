"""Base hardware orchestration for the master-satellite architecture."""

from artiq.experiment import kernel, delay, ms, us

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

    def __init__(
        self,
        experiment,
        experiment_mode,
        which_node=None,
        satellite_wait_attempts=100,
        satellite_ready_poll_interval=100 * ms,
    ):
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

        if satellite_wait_attempts < 1:
            raise ValueError("satellite_wait_attempts must be positive.")

        self.experiment = experiment
        self.experiment_mode = experiment_mode
        self.which_node = which_node
        self.active_nodes = active_nodes
        self.node2_active = "Node2" in active_nodes
        self.satellite_wait_attempts = int(satellite_wait_attempts)
        self.satellite_ready_poll_interval = satellite_ready_poll_interval

        self.node_resolvers = {}
        self.shared_spcm_resolver = None
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

        for device_name in ("core", "core_dma", "scheduler"):
            self.experiment.setattr_device(device_name)

        for node in self.active_nodes:
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
        self._built = True

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

        if self.experiment_mode == "single_node" and self.which_node == "Node2":
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
        if self._prepared:
            raise RuntimeError("BaseExperimentMasterSatellite.prepare() called twice.")

        for node in self.active_nodes:
            resolver = self.node_resolvers[node]
            self._publish_compatibility_default_attributes(node, resolver)
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

    def _publish_compatibility_default_attributes(self, node, resolver):
        if self.experiment_mode != "single_node":
            return

        default_names = set()
        for defaults in resolver.dds_defaults.values():
            default_names.add(defaults["frequency"])
            default_names.add(defaults["power"])

        for compatibility_name in default_names:
            authoritative_name = self._default_attribute_name(
                node, compatibility_name
            )
            try:
                authoritative_value = getattr(
                    self.experiment, authoritative_name
                )
            except AttributeError as error:
                raise RuntimeError(
                    f"Missing authoritative master-satellite variable "
                    f"{authoritative_name!r} for {node}."
                ) from error
            setattr(self.experiment, compatibility_name, authoritative_value)

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
        self, turn_off_dds_channels=True, turn_off_zotinos=True
    ):
        """Initialize active hardware using one core reset and ordered access."""
        if not self._prepared:
            raise RuntimeError("prepare() must be called before initialize_hardware().")

        # This is the only core reset owned by the master-satellite base.
        self.experiment.core.reset()

        if self.node2_active:
            self._wait_for_satellite()
        else:
            self.experiment.core.break_realtime()

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
            devices[i].set_att(0.0)
            amplitude = (2 * 50 * 10 ** (powers[i] / 10 - 3)) ** 0.5
            devices[i].set(
                frequency=frequencies[i], amplitude=amplitude
            )
            if turn_off_outputs:
                devices[i].sw.off()
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
