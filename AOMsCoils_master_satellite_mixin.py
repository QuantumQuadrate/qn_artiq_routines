"""Shared implementation for node-local master-satellite manual controls."""

from artiq.experiment import BooleanValue, NumberValue, delay, kernel, ms

from subroutines.k10cr1_functions import (
    go_to_home,
    move_by_deg,
    move_to_target_deg,
)
from utilities.BaseExperiment_master_satellite import BaseExperimentMasterSatellite
from utilities.conversions import dB_to_V_kernel as dB_to_V


class _AOMsCoilsMasterSatelliteMixin:
    """Common implementation; public node experiments live in separate files.

    This class deliberately does not inherit EnvExperiment and therefore must
    not appear as a separate ARTIQ Explorer experiment.
    """

    NODE = None
    COIL_CHANNELS = ()

    K10CR1_CONTROL_ARGUMENTS = (
        "go_to_home_780HWP",
        "go_to_home_780QWP",
        "go_to_target_780HWP",
        "go_to_target_780QWP",
        "move_780HWP_by",
        "move_780QWP_by",
        "go_to_home_852HWP",
        "go_to_home_852QWP",
        "go_to_target_852HWP",
        "go_to_target_852QWP",
        "move_852HWP_by",
        "move_852QWP_by",
        "go_to_optimized_852_settings",
    )

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
        # Coils follow the standalone AOMsCoils semantics: unless disabled,
        # all four coils are driven to the node's persistent MOT calibration
        # voltages. There are no manual voltage entries.
        self.setattr_argument(
            "disable_coils", BooleanValue(False), f"{self.NODE} coils"
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
        for control_name in self.K10CR1_CONTROL_ARGUMENTS:
            group = (
                "K10CR1 780 waveplates"
                if "780" in control_name
                else "K10CR1 852 waveplates"
            )
            self.setattr_argument(control_name, BooleanValue(False), group)

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

        # MOT coil voltages come from the node's persistent MOT calibration
        # through the single-node projection; disable_coils selects 0 V.
        self._coil_voltages = [
            self.AZ_bottom_volts_MOT,
            self.AZ_top_volts_MOT,
            self.AX_volts_MOT,
            self.AY_volts_MOT,
        ]

        # One controller process owns all eight rotators; the node identity
        # lives in the axis nickname. The NDSP client opens its connection
        # the moment the device is requested, so it is bound only when a
        # waveplate action is actually selected.
        self._k10cr1_requested = any(
            bool(getattr(self, name)) for name in self.K10CR1_CONTROL_ARGUMENTS
        )
        if self._k10cr1_requested:
            self.setattr_device("k10cr1_ndsp")
        self._axis_780_HWP = f"780_HWP_{self.NODE}"
        self._axis_780_QWP = f"780_QWP_{self.NODE}"
        self._axis_852_HWP = f"852_HWP_{self.NODE}"
        self._axis_852_QWP = f"852_QWP_{self.NODE}"

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
        if self.disable_coils:
            self.zotino0.set_dac(
                [0.0, 0.0, 0.0, 0.0], channels=self.COIL_CHANNELS
            )
        else:
            self.zotino0.set_dac(
                self._coil_voltages, channels=self.COIL_CHANNELS
            )
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
    def k10cr1_operations(self):
        """Rotate only this node's waveplates; axis names carry the node.

        Unlike the standalone utility, the 852 target moves are gated on the
        852 booleans (the original gated them on the 780 booleans).
        """
        self.core.reset()
        delay(10 * ms)
        if self.go_to_home_780HWP:
            go_to_home(self, self._axis_780_HWP)
        if self.go_to_home_780QWP:
            go_to_home(self, self._axis_780_QWP)

        if self.go_to_target_780HWP:
            move_to_target_deg(
                self, name=self._axis_780_HWP, target_deg=self.target_780_HWP
            )
        if self.go_to_target_780QWP:
            move_to_target_deg(
                self, name=self._axis_780_QWP, target_deg=self.target_780_QWP
            )

        if self.move_780HWP_by:
            move_by_deg(
                self, name=self._axis_780_HWP, target_deg=self.move_780_HWP_by
            )
        if self.move_780QWP_by:
            move_by_deg(
                self, name=self._axis_780_QWP, target_deg=self.move_780_QWP_by
            )

        if self.go_to_home_852HWP:
            go_to_home(self, self._axis_852_HWP)
        if self.go_to_home_852QWP:
            go_to_home(self, self._axis_852_QWP)

        if self.go_to_target_852HWP:
            move_to_target_deg(
                self, name=self._axis_852_HWP, target_deg=self.target_852_HWP
            )
        if self.go_to_target_852QWP:
            move_to_target_deg(
                self, name=self._axis_852_QWP, target_deg=self.target_852_QWP
            )

        if self.move_852HWP_by:
            move_by_deg(
                self, name=self._axis_852_HWP, target_deg=self.move_852_HWP_by
            )
        if self.move_852QWP_by:
            move_by_deg(
                self, name=self._axis_852_QWP, target_deg=self.move_852_QWP_by
            )

        if self.go_to_optimized_852_settings:
            move_to_target_deg(
                self,
                name=self._axis_852_HWP,
                target_deg=self.best_852HWP_to_max,
            )
            move_to_target_deg(
                self,
                name=self._axis_852_QWP,
                target_deg=self.best_852QWP_to_max,
            )

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
        if self._k10cr1_requested:
            self.k10cr1_operations()
