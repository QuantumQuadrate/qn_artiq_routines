"""Shared implementation for node-local master-satellite manual controls."""

from artiq.experiment import BooleanValue, NumberValue, delay, kernel, ms

from subroutines.k10cr1_functions import (
    go_to_home,
    move_by_deg,
    move_to_target_deg,
)
from utilities.BaseExperiment_master_satellite import (
    BaseExperimentMasterSatellite,
    _DatasetRedirectMixin,
)
from utilities.conversions import dB_to_V_kernel as dB_to_V


class _AOMsCoilsMasterSatelliteMixin(_DatasetRedirectMixin):
    """Common implementation; public node experiments live in separate files.

    This class deliberately does not inherit EnvExperiment and therefore must
    not appear as a separate ARTIQ Explorer experiment.
    """

    NODE = None
    COIL_CHANNELS = ()
    # Tests inject a fake stabilizer class here; None selects the real
    # AOMPowerStabilizer from subroutines/aom_feedback.py.
    _stabilizer_factory = None

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
        if self.NODE == "Node2":
            # Standalone bob-only combined modes: DDS and optical gate driven
            # together. Applied last, so on Node2 they are authoritative over
            # the individual GRIN controls (their OFF branches override).
            self.setattr_argument(
                "Node2_GRIN1_AOM_ON",
                BooleanValue(False),
                "Node2 GRIN1/GRIN2 combined modes",
            )
            self.setattr_argument(
                "Node2_GRIN2_AOM_ON",
                BooleanValue(False),
                "Node2 GRIN1/GRIN2 combined modes",
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
        self.setattr_argument(
            "run_laser_feedback",
            BooleanValue(False),
            "Laser power stabilization",
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

        # The stabilizer reads its per-node feedback config via the legacy
        # node name and persists through the node-suffixed datasets.
        self.which_node = self.base.NODE_LEGACY_NAMES[self.NODE]
        self.base.prepare_laser_stabilizer(
            stabilizer_factory=self._stabilizer_factory
        )
        self._all_stabilized_aoms_requested = bool(
            self.AOM_A1_ON and self.AOM_A2_ON and self.AOM_A3_ON
            and self.AOM_A4_ON and self.AOM_A5_ON and self.AOM_A6_ON
            and self.cooling_DP_ON
        )

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

        # f/p_GRIN2_excitation project on both nodes, so the combined-mode
        # kernel compiles everywhere; it only runs on Node2.
        self._apply_node2_combined_modes = self.NODE == "Node2"
        self._node2_grin1_combined = bool(
            getattr(self, "Node2_GRIN1_AOM_ON", False)
        )
        self._node2_grin2_combined = bool(
            getattr(self, "Node2_GRIN2_AOM_ON", False)
        )
        self._grin2_excitation_frequency = self.f_GRIN2_excitation
        self._grin2_excitation_amplitude_dB = self.p_GRIN2_excitation

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
    def _apply_node2_combined_grin_state(self):
        """Replicate the standalone bob-only combined GRIN modes.

        Applied after the individual DDS/gate controls, exactly as in the
        standalone AOMsCoils ordering, so on Node2 these switches are
        authoritative: their OFF branches override GRIN1and2_ON, the GRIN
        gate controls, and the D1 DDS (D1_pumping_DP_ON on Node2 therefore
        effectively drives only ttl_D1_pumping, as in the original).
        """
        if self._node2_grin1_combined:
            self.GRIN1and2_dds.sw.on()
            self.ttl_GRIN1_switch.off()
        else:
            self.GRIN1and2_dds.sw.off()
            self.ttl_GRIN1_switch.on()
        delay(1 * ms)
        if self._node2_grin2_combined:
            self.dds_D1_pumping_DP.set(
                frequency=self._grin2_excitation_frequency,
                amplitude=dB_to_V(self._grin2_excitation_amplitude_dB),
            )
            self.dds_D1_pumping_DP.sw.on()
            self.ttl_GRIN2_switch.off()
        else:
            self.dds_D1_pumping_DP.sw.off()
            self.ttl_GRIN2_switch.on()
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
    def run_feedback(self):
        """Feed back (or monitor) when every stabilized AOM is requested on.

        Mirrors the standalone utility: the stabilizer only runs when all six
        fiber AOMs and the cooling DP are requested on, and the requested DDS
        and optical-gate states are re-applied afterward because the
        stabilizer switches AOMs during its measurement phase.
        """
        self.core.reset()
        if self._all_stabilized_aoms_requested:
            if self.run_laser_feedback:
                self.laser_stabilizer.run()
                delay(1 * ms)
                self._apply_dds_state()
                delay(1 * ms)
                self._apply_active_low_ttl_state()
            else:
                delay(10 * ms)
                self.laser_stabilizer.monitor()
                delay(1 * ms)
                self._apply_dds_state()
                delay(1 * ms)
                self._apply_active_low_ttl_state()

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
        if self._apply_node2_combined_modes:
            self._apply_node2_combined_grin_state()
        self._apply_microwave_state()
        self._apply_coil_state()

    def run(self):
        self.apply_manual_state()
        self.run_feedback()
        if self._k10cr1_requested:
            self.k10cr1_operations()
