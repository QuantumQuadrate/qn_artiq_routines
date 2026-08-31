"""Shared implementation for node-local master-satellite manual controls."""

from artiq.experiment import BooleanValue, NumberValue, delay, kernel, ms

from utilities.BaseExperiment_master_satellite import BaseExperimentMasterSatellite
from utilities.conversions import dB_to_V_kernel as dB_to_V


class _AOMsCoilsMasterSatelliteMixin:
    """Common implementation; public node experiments live in separate files.

    This class deliberately does not inherit EnvExperiment and therefore must
    not appear as a separate ARTIQ Explorer experiment.
    """

    NODE = None
    COIL_CHANNELS = ()

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

    def _build_node_manual_controls(self):
        if self.NODE not in ("Node1", "Node2") or len(self.COIL_CHANNELS) != 4:
            raise RuntimeError(
                "Public AOMsCoils experiment must define NODE and four "
                "explicit COIL_CHANNELS."
            )

        # Examination-safe superset registration is narrowed to this node's
        # hardware lifecycle during prepare(). The other node is never
        # initialized or assigned an output state by this experiment.
        self.base = BaseExperimentMasterSatellite(experiment=self)
        self.base.build()

        for control_name, _ in self.DDS_CONTROLS:
            self.setattr_argument(
                control_name, BooleanValue(False), f"{self.NODE} DDS/AOM outputs"
            )
        for control_name, _ in self.ACTIVE_LOW_TTL_CONTROLS:
            self.setattr_argument(
                control_name,
                BooleanValue(False),
                f"{self.NODE} TTL-controlled AOM gates",
            )
        for coil_name in self.COIL_NAMES:
            self.setattr_argument(
                f"{coil_name}_coil_ON",
                BooleanValue(False),
                f"{self.NODE} coils",
            )
            self.setattr_argument(
                f"{coil_name}_voltage",
                NumberValue(
                    0.0, unit="V", min=-10.0, max=10.0, ndecimals=3
                ),
                f"{self.NODE} run-local coil voltage (safe 0 V default)",
            )
        self.setattr_argument(
            "microwave_dds_ON",
            BooleanValue(False),
            f"{self.NODE} microwave/RF outputs",
        )
        self.setattr_argument(
            "MW_RF_dds_ON",
            BooleanValue(False),
            f"{self.NODE} microwave/RF outputs",
        )
        self.setattr_argument(
            "yes_Im_sure_I_want_MW_or_RF_dds_ON",
            BooleanValue(False),
            "Microwave/RF safety confirmation",
        )

    def prepare(self):
        self.base.configure_execution("single_node", self.NODE)
        self.base.prepare()

        resolver = self.base.node_resolvers[self.NODE]
        self._manual_dds = []
        self._manual_dds_enabled = []
        self._manual_dds_frequencies = []
        self._manual_dds_powers = []
        for control_name, logical_alias in self.DDS_CONTROLS:
            binding = next(
                item
                for item in resolver.dds_bindings
                if item["logical_alias"] == logical_alias
            )
            self._manual_dds.append(
                getattr(self, binding["experiment_attribute"])
            )
            self._manual_dds_enabled.append(getattr(self, control_name))
            self._manual_dds_frequencies.append(
                getattr(self, binding["frequency_attribute"])
            )
            self._manual_dds_powers.append(
                getattr(self, binding["power_attribute"])
            )

        self._manual_ttls = []
        self._manual_ttl_enabled = []
        for control_name, logical_alias in self.ACTIVE_LOW_TTL_CONTROLS:
            self._manual_ttls.append(getattr(self, logical_alias))
            enabled = getattr(self, control_name)
            if self.NODE == "Node1" and logical_alias == "ttl_GRIN1_switch":
                enabled = enabled or self.D1_pumping_DP_ON
            self._manual_ttl_enabled.append(enabled)
        if self.NODE == "Node2":
            self._manual_ttls.append(self.ttl_D1_pumping)
            self._manual_ttl_enabled.append(self.D1_pumping_DP_ON)

        self._coil_enabled = [
            getattr(self, f"{name}_coil_ON") for name in self.COIL_NAMES
        ]
        self._coil_voltages = [
            getattr(self, f"{name}_voltage") for name in self.COIL_NAMES
        ]

    @kernel
    def _force_selected_node_off(self):
        for dds in self._manual_dds:
            dds.sw.off()
            delay(1 * ms)
        for ttl in self._manual_ttls:
            ttl.on()
            delay(1 * ms)
        self.zotino0.set_dac(
            [0.0, 0.0, 0.0, 0.0], channels=self.COIL_CHANNELS
        )
        delay(1 * ms)

    @kernel
    def _apply_dds_state(self):
        for i in range(len(self._manual_dds)):
            if self._manual_dds_enabled[i]:
                self._manual_dds[i].set(
                    frequency=self._manual_dds_frequencies[i],
                    amplitude=dB_to_V(self._manual_dds_powers[i]),
                )
                self._manual_dds[i].sw.on()
            else:
                self._manual_dds[i].sw.off()
            delay(1 * ms)

    @kernel
    def _apply_active_low_ttl_state(self):
        for i in range(len(self._manual_ttls)):
            if self._manual_ttl_enabled[i]:
                self._manual_ttls[i].off()
            else:
                self._manual_ttls[i].on()
            delay(1 * ms)

    @kernel
    def _apply_coil_state(self):
        requested = [0.0, 0.0, 0.0, 0.0]
        for i in range(4):
            if self._coil_enabled[i]:
                requested[i] = self._coil_voltages[i]
        self.zotino0.set_dac(requested, channels=self.COIL_CHANNELS)
        delay(1 * ms)

    @kernel
    def _apply_microwave_state(self):
        allow = self.yes_Im_sure_I_want_MW_or_RF_dds_ON
        if allow and self.microwave_dds_ON:
            self.dds_microwaves.sw.on()
            self.ttl_microwave_switch.off()
        else:
            self.dds_microwaves.sw.off()
            self.ttl_microwave_switch.on()
        if allow and self.MW_RF_dds_ON:
            self.dds_MW_RF.sw.on()
        else:
            self.dds_MW_RF.sw.off()
        delay(1 * ms)

    @kernel
    def apply_manual_state(self):
        """Initialize and modify only this experiment's selected node."""
        self.base.initialize_hardware(
            turn_off_dds_channels=True,
            turn_off_zotinos=True,
            initialize_spcms=False,
            reset_core=False,
        )
        self._force_selected_node_off()
        self._apply_dds_state()
        self._apply_active_low_ttl_state()
        self._apply_microwave_state()
        self._apply_coil_state()

    def run(self):
        self.apply_manual_state()
