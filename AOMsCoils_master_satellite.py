"""Deterministic manual AOM/DDS and coil control for master-satellite hardware.

This utility intentionally binds both physical nodes at all times.  Its
``which_node`` argument is a desired-hardware-state gate, not a device-binding
selection: outputs belonging to an unselected node are actively disabled.
"""

from artiq.experiment import (
    BooleanValue,
    EnumerationValue,
    EnvExperiment,
    NumberValue,
    delay,
    kernel,
    ms,
)

from utilities.BaseExperiment_master_satellite import (
    BaseExperimentMasterSatellite,
)
from utilities.conversions import dB_to_V_kernel as dB_to_V


class AOMsCoils_master_satellite(EnvExperiment):
    """AOMsCoils_master_satellite

    Manually apply an explicit, safely gated state to both nodes.
    """

    VALID_NODE_SELECTIONS = ("node1", "node2", "two_nodes")

    DDS_CONTROLS = (
        ("FORT_AOM_ON", "dds_FORT"),
        ("cooling_DP_ON", "dds_cooling_DP"),
        ("AOM_A1_ON", "dds_AOM_A1"),
        ("AOM_A2_ON", "dds_AOM_A2"),
        ("AOM_A3_ON", "dds_AOM_A3"),
        ("AOM_A4_ON", "dds_AOM_A4"),
        ("AOM_A5_ON", "dds_AOM_A5"),
        ("AOM_A6_ON", "dds_AOM_A6"),
        ("D1_pumping_DP_ON", "dds_D1_pumping_DP"),
        ("GRIN1and2_ON", "GRIN1and2_dds"),
    )

    ACTIVE_LOW_TTL_CONTROLS = (
        ("Repump_AOM_switch_ON", "ttl_repump_switch"),
        ("pumping_repump_switch_ON", "ttl_pumping_repump_switch"),
        ("Excitation0_AOM_switch_ON", "ttl_exc0_switch"),
        ("GRIN1_AOM_switch_ON", "ttl_GRIN1_switch"),
        ("GRIN2_AOM_switch_ON", "ttl_GRIN2_switch"),
    )

    COIL_NAMES = ("AZ_bottom", "AZ_top", "AX", "AY")
    COIL_CHANNELS_NODE1 = (0, 1, 13, 14)
    COIL_CHANNELS_NODE2 = (0, 1, 2, 3)

    def build(self):
        # Both nodes must remain accessible even in node1/node2 selections so
        # that the unselected node can be actively forced to a safe state.
        # Execution configuration is deferred so repository examination does
        # not require persistent calibration datasets to exist.
        self.base = BaseExperimentMasterSatellite(experiment=self)
        self.base.build()

        self.setattr_argument(
            "which_node",
            EnumerationValue(self.VALID_NODE_SELECTIONS, default="node1"),
            "Desired hardware state",
        )

        for node in ("Node1", "Node2"):
            for control_name, _ in self.DDS_CONTROLS:
                self.setattr_argument(
                    f"{control_name}_{node}",
                    BooleanValue(False),
                    f"{node} DDS/AOM outputs",
                )
            for control_name, _ in self.ACTIVE_LOW_TTL_CONTROLS:
                self.setattr_argument(
                    f"{control_name}_{node}",
                    BooleanValue(False),
                    f"{node} TTL-controlled AOM gates",
                )
            for coil_name in self.COIL_NAMES:
                self.setattr_argument(
                    f"{coil_name}_coil_ON_{node}",
                    BooleanValue(False),
                    f"{node} coils",
                )
                self.setattr_argument(
                    f"{coil_name}_voltage_{node}",
                    NumberValue(
                        0.0,
                        unit="V",
                        min=-10.0,
                        max=10.0,
                        ndecimals=3,
                    ),
                    f"{node} run-local coil voltages (safe 0 V default)",
                )

            self.setattr_argument(
                f"microwave_dds_ON_{node}",
                BooleanValue(False),
                f"{node} microwave/RF outputs",
            )
            self.setattr_argument(
                f"MW_RF_dds_ON_{node}",
                BooleanValue(False),
                f"{node} microwave/RF outputs",
            )

        self.setattr_argument(
            "yes_Im_sure_I_want_MW_or_RF_dds_ON",
            BooleanValue(False),
            "Microwave/RF safety confirmation",
        )

    def prepare(self):
        if self.which_node not in self.VALID_NODE_SELECTIONS:
            raise ValueError(
                f"Unsupported which_node {self.which_node!r}; expected "
                "'node1', 'node2', or 'two_nodes'."
            )

        # This utility always needs both nodes physically available even when
        # one node is selected, because the other node is actively forced OFF.
        # Strict authoritative dataset loading happens here, before hardware.
        self.base.configure_execution("two_nodes")
        self.base.prepare()

        self.node1_active = self.which_node in ("node1", "two_nodes")
        self.node2_active = self.which_node in ("node2", "two_nodes")

        self._prepare_node_state("Node1")
        self._prepare_node_state("Node2")

    def _prepare_node_state(self, node):
        """Capture run-local GUI state and authoritative DDS defaults."""
        resolver = self.base.node_resolvers[node]
        suffix = node.lower()

        dds_devices = []
        dds_enabled = []
        dds_frequencies = []
        dds_powers = []
        for control_name, logical_alias in self.DDS_CONTROLS:
            binding = next(
                item
                for item in resolver.dds_bindings
                if item["logical_alias"] == logical_alias
            )
            dds_devices.append(
                getattr(self, binding["experiment_attribute"])
            )
            dds_enabled.append(getattr(self, f"{control_name}_{node}"))
            dds_frequencies.append(
                getattr(self, binding["frequency_attribute"])
            )
            dds_powers.append(getattr(self, binding["power_attribute"]))

        setattr(self, f"_manual_dds_{suffix}", dds_devices)
        setattr(self, f"_manual_dds_enabled_{suffix}", dds_enabled)
        setattr(self, f"_manual_dds_frequencies_{suffix}", dds_frequencies)
        setattr(self, f"_manual_dds_powers_{suffix}", dds_powers)

        ttl_devices = []
        ttl_enabled = []
        for control_name, logical_alias in self.ACTIVE_LOW_TTL_CONTROLS:
            ttl_devices.append(getattr(self, f"{logical_alias}_{node}"))
            enabled = getattr(self, f"{control_name}_{node}")
            # Node1 D1 pumping shares the GRIN1 active-low optical gate.
            if node == "Node1" and logical_alias == "ttl_GRIN1_switch":
                enabled = enabled or self.D1_pumping_DP_ON_Node1
            ttl_enabled.append(enabled)
        if node == "Node2":
            # Bob's standalone D1 optical path has a dedicated active-low
            # gate.  It follows the D1 DDS request and has no second GUI state.
            ttl_devices.append(self.ttl_D1_pumping_Node2)
            ttl_enabled.append(self.D1_pumping_DP_ON_Node2)
        setattr(self, f"_manual_ttls_{suffix}", ttl_devices)
        setattr(self, f"_manual_ttl_enabled_{suffix}", ttl_enabled)

        coil_enabled = [
            getattr(self, f"{name}_coil_ON_{node}")
            for name in self.COIL_NAMES
        ]
        coil_voltages = [
            getattr(self, f"{name}_voltage_{node}")
            for name in self.COIL_NAMES
        ]
        setattr(self, f"_coil_enabled_{suffix}", coil_enabled)
        setattr(self, f"_coil_voltages_{suffix}", coil_voltages)

    @kernel
    def _force_node_off(self, dds_devices, ttl_devices, zotino, channels):
        for dds in dds_devices:
            dds.sw.off()
            delay(1 * ms)
        for ttl in ttl_devices:
            # These optical gates are active-low; high is the safe OFF state.
            ttl.on()
            delay(1 * ms)
        zotino.set_dac([0.0, 0.0, 0.0, 0.0], channels=channels)
        delay(1 * ms)

    @kernel
    def _apply_dds_state(
        self, devices, enabled, frequencies, powers, node_active
    ):
        for i in range(len(devices)):
            if node_active and enabled[i]:
                devices[i].set(
                    frequency=frequencies[i], amplitude=dB_to_V(powers[i])
                )
                devices[i].sw.on()
            else:
                devices[i].sw.off()
            delay(1 * ms)

    @kernel
    def _apply_active_low_ttl_state(self, devices, enabled, node_active):
        for i in range(len(devices)):
            if node_active and enabled[i]:
                devices[i].off()
            else:
                devices[i].on()
            delay(1 * ms)

    @kernel
    def _apply_coil_state(
        self, zotino, channels, enabled, voltages, node_active
    ):
        requested = [0.0, 0.0, 0.0, 0.0]
        for i in range(4):
            if node_active and enabled[i]:
                requested[i] = voltages[i]
        zotino.set_dac(requested, channels=channels)
        delay(1 * ms)

    @kernel
    def _apply_microwave_state(self, node):
        # Microwave paths require both an explicit per-node request and the
        # run-wide safety confirmation.  The TTL switch is active-low.
        if node == 1:
            active = self.node1_active
            microwave_on = self.microwave_dds_ON_Node1
            rf_on = self.MW_RF_dds_ON_Node1
            microwave = self.dds_microwaves_Node1
            rf = self.dds_MW_RF_Node1
            switch = self.ttl_microwave_switch_Node1
        else:
            active = self.node2_active
            microwave_on = self.microwave_dds_ON_Node2
            rf_on = self.MW_RF_dds_ON_Node2
            microwave = self.dds_microwaves_Node2
            rf = self.dds_MW_RF_Node2
            switch = self.ttl_microwave_switch_Node2

        allow = active and self.yes_Im_sure_I_want_MW_or_RF_dds_ON
        if allow and microwave_on:
            microwave.sw.on()
            switch.off()
        else:
            microwave.sw.off()
            switch.on()
        if allow and rf_on:
            rf.sw.on()
        else:
            rf.sw.off()
        delay(1 * ms)

    @kernel
    def apply_manual_state(self):
        """Initialize safely, force both nodes OFF, then apply GUI state."""
        self.base.initialize_hardware(
            turn_off_dds_channels=True, turn_off_zotinos=True
        )

        # This explicit second OFF pass is intentional.  It documents and
        # enforces the distinction between raw TTL wiring and desired optical
        # state, and guarantees hard gating even if Base initialization evolves.
        self._force_node_off(
            self._manual_dds_node1,
            self._manual_ttls_node1,
            self.zotino0_Node1,
            self.COIL_CHANNELS_NODE1,
        )
        self._force_node_off(
            self._manual_dds_node2,
            self._manual_ttls_node2,
            self.zotino0_Node2,
            self.COIL_CHANNELS_NODE2,
        )

        self._apply_dds_state(
            self._manual_dds_node1,
            self._manual_dds_enabled_node1,
            self._manual_dds_frequencies_node1,
            self._manual_dds_powers_node1,
            self.node1_active,
        )
        self._apply_dds_state(
            self._manual_dds_node2,
            self._manual_dds_enabled_node2,
            self._manual_dds_frequencies_node2,
            self._manual_dds_powers_node2,
            self.node2_active,
        )
        self._apply_active_low_ttl_state(
            self._manual_ttls_node1,
            self._manual_ttl_enabled_node1,
            self.node1_active,
        )
        self._apply_active_low_ttl_state(
            self._manual_ttls_node2,
            self._manual_ttl_enabled_node2,
            self.node2_active,
        )
        self._apply_microwave_state(1)
        self._apply_microwave_state(2)
        self._apply_coil_state(
            self.zotino0_Node1,
            self.COIL_CHANNELS_NODE1,
            self._coil_enabled_node1,
            self._coil_voltages_node1,
            self.node1_active,
        )
        self._apply_coil_state(
            self.zotino0_Node2,
            self.COIL_CHANNELS_NODE2,
            self._coil_enabled_node2,
            self._coil_voltages_node2,
            self.node2_active,
        )

    def run(self):
        self.apply_manual_state()
