"""AOMsCoils_master_satellite_Node2

Manual control of Node2's AOMs, optical gates, microwaves, MOT coils and
waveplates, without changing anything on Node1. The master-satellite
counterpart of standalone/AOMsCoils.py on bob, and written the same way:
one block per GUI box.

This file is complete on its own. AOMsCoils_master_satellite_Node1.py is a
separate file for Node1; the two share no code, so a change for one node's
hardware cannot reach the other.

The TTL gates are ACTIVE LOW: off() lets the light (or RF) through, on()
blocks it.
"""

from artiq.experiment import BooleanValue, EnvExperiment, delay, kernel, ms

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


class AOMsCoils_master_satellite_Node2(_DatasetRedirectMixin, EnvExperiment):
    """AOMsCoils_master_satellite_Node2

    Manually control only Node2 hardware without changing Node1 outputs.
    """

    NODE = "Node2"
    # Tests inject a fake stabilizer class here; None selects the real
    # AOMPowerStabilizer from subroutines/aom_feedback.py.
    _stabilizer_factory = None

    def build(self):
        # build() registers both nodes' devices for the repository scan;
        # prepare() narrows the run to Node2. Node1 is never initialized or
        # given an output state by this experiment.
        self.base = BaseExperimentMasterSatellite(experiment=self)
        self.base.build()

        group = "Node2 DDS/AOM outputs"
        self.setattr_argument("FORT_AOM_ON", BooleanValue(False), group)
        self.setattr_argument("cooling_DP_ON", BooleanValue(False), group)
        self.setattr_argument("AOM_A1_ON", BooleanValue(False), group)
        self.setattr_argument("AOM_A2_ON", BooleanValue(False), group)
        self.setattr_argument("AOM_A3_ON", BooleanValue(False), group)
        self.setattr_argument("AOM_A4_ON", BooleanValue(False), group)
        self.setattr_argument("AOM_A5_ON", BooleanValue(False), group)
        self.setattr_argument("AOM_A6_ON", BooleanValue(False), group)
        self.setattr_argument("GRIN1and2_ON", BooleanValue(False), group)

        group = "Node2 TTL-controlled AOM gates"
        self.setattr_argument("D1_pumping_switch_ON", BooleanValue(False), group)
        self.setattr_argument("Repump_AOM_switch_ON", BooleanValue(False), group)
        self.setattr_argument("pumping_repump_switch_ON", BooleanValue(False), group)
        self.setattr_argument("Excitation0_AOM_switch_ON", BooleanValue(False), group)
        self.setattr_argument("GRIN1_AOM_switch_ON", BooleanValue(False), group)
        self.setattr_argument("GRIN2_AOM_switch_ON", BooleanValue(False), group)

        # Node2 only, as in standalone bob: each box drives a GRIN DDS and
        # its gate together, and has the final say over them (turn_on_AOMs).
        group = "Node2 GRIN1/GRIN2 combined modes"
        self.setattr_argument("Node2_GRIN1_AOM_ON", BooleanValue(False), group)
        self.setattr_argument("Node2_GRIN2_AOM_ON", BooleanValue(False), group)

        # Unless disabled, the four coils go to Node2's persistent MOT
        # calibration voltages. There are no manual voltage entries.
        self.setattr_argument("disable_coils", BooleanValue(False), "Node2 coils")

        group = "Node2 microwave/RF outputs"
        self.setattr_argument("microwave_dds_ON", BooleanValue(False), group)
        self.setattr_argument("MW_RF_dds_ON", BooleanValue(False), group)
        self.setattr_argument(
            "yes_Im_sure_I_want_MW_or_RF_dds_ON",
            BooleanValue(False),
            "Microwave/RF safety confirmation",
        )

        self.setattr_argument(
            "run_laser_feedback", BooleanValue(False), "Laser power stabilization"
        )

        group = "K10CR1 780 waveplates"
        self.setattr_argument("go_to_home_780HWP", BooleanValue(False), group)
        self.setattr_argument("go_to_home_780QWP", BooleanValue(False), group)
        self.setattr_argument("go_to_target_780HWP", BooleanValue(False), group)
        self.setattr_argument("go_to_target_780QWP", BooleanValue(False), group)
        self.setattr_argument("move_780HWP_by", BooleanValue(False), group)
        self.setattr_argument("move_780QWP_by", BooleanValue(False), group)

        group = "K10CR1 852 waveplates"
        self.setattr_argument("go_to_home_852HWP", BooleanValue(False), group)
        self.setattr_argument("go_to_home_852QWP", BooleanValue(False), group)
        self.setattr_argument("go_to_target_852HWP", BooleanValue(False), group)
        self.setattr_argument("go_to_target_852QWP", BooleanValue(False), group)
        self.setattr_argument("move_852HWP_by", BooleanValue(False), group)
        self.setattr_argument("move_852QWP_by", BooleanValue(False), group)
        self.setattr_argument(
            "go_to_optimized_852_settings", BooleanValue(False), group
        )

    def prepare(self):
        # Load only Node2's variables (plus the globals). The bare device
        # names used below (dds_FORT, ttl_repump_switch, zotino0, ...) are
        # Node2's from here on.
        self.base.configure_execution("single_node", self.NODE)
        self.base.prepare()

        # The stabilizer reads its config through the legacy node name
        # ("bob") and persists through the _Node2 datasets.
        self.which_node = self.base.NODE_LEGACY_NAMES[self.NODE]
        self.base.prepare_laser_stabilizer(
            stabilizer_factory=self._stabilizer_factory
        )

        # As in standalone, feedback runs only when the six fiber AOMs and
        # the cooling DP are all requested on.
        self._all_stabilized_aoms_requested = bool(
            self.AOM_A1_ON and self.AOM_A2_ON and self.AOM_A3_ON
            and self.AOM_A4_ON and self.AOM_A5_ON and self.AOM_A6_ON
            and self.cooling_DP_ON
        )

        self._k10cr1_requested = bool(
            self.go_to_home_780HWP or self.go_to_home_780QWP
            or self.go_to_target_780HWP or self.go_to_target_780QWP
            or self.move_780HWP_by or self.move_780QWP_by
            or self.go_to_home_852HWP or self.go_to_home_852QWP
            or self.go_to_target_852HWP or self.go_to_target_852QWP
            or self.move_852HWP_by or self.move_852QWP_by
            or self.go_to_optimized_852_settings
        )

    @kernel
    def everything_off(self):
        """The AOM DDSes, optical gates and coils: RF off, gates blocked,
        coils at 0 V. initialize_hardware() has already done this (and turned
        the microwaves off); this pass makes the starting point explicit
        before the boxes are applied.
        """
        self.dds_FORT.sw.off()
        delay(1 * ms)
        self.dds_cooling_DP.sw.off()
        delay(1 * ms)
        self.dds_AOM_A1.sw.off()
        delay(1 * ms)
        self.dds_AOM_A2.sw.off()
        delay(1 * ms)
        self.dds_AOM_A3.sw.off()
        delay(1 * ms)
        self.dds_AOM_A4.sw.off()
        delay(1 * ms)
        self.dds_AOM_A5.sw.off()
        delay(1 * ms)
        self.dds_AOM_A6.sw.off()
        delay(1 * ms)
        self.dds_GRIN2.sw.off()
        delay(1 * ms)
        self.dds_GRIN1.sw.off()
        delay(1 * ms)

        self.ttl_repump_switch.on()
        delay(1 * ms)
        self.ttl_pumping_repump_switch.on()
        delay(1 * ms)
        self.ttl_exc0_switch.on()
        delay(1 * ms)
        self.ttl_GRIN1_switch.on()
        delay(1 * ms)
        self.ttl_GRIN2_switch.on()
        delay(1 * ms)
        self.ttl_D1_pumping_switch.on()
        delay(1 * ms)

        self.zotino0.set_dac([0.0, 0.0, 0.0, 0.0], channels=self.coil_channels)
        delay(1 * ms)

    @kernel
    def turn_on_AOMs(self):
        """Turn on what the boxes ask for and turn off the rest.

        Run after everything_off(), and again after laser feedback, which
        switches AOMs during its measurements. As in standalone, this only
        switches, except that the combined GRIN2 mode sets GRIN2's
        excitation tone. initialize_hardware() has set every DDS to Node2's
        frequency and power, and laser feedback corrects the power of each
        AOM it stabilizes (fast_feedback_dds_list) in the DDS itself, so
        after feedback the corrected power is what comes back on.

        The two combined GRIN modes come last, as in standalone bob, and
        decide how both GRIN DDSes and both GRIN gates end up.
        """
        ### DDS / AOM outputs: switched only (see above)
        if self.FORT_AOM_ON:
            self.dds_FORT.sw.on()
        else:
            self.dds_FORT.sw.off()
        delay(1 * ms)

        if self.cooling_DP_ON:
            self.dds_cooling_DP.sw.on()
        else:
            self.dds_cooling_DP.sw.off()
        delay(1 * ms)

        if self.AOM_A1_ON:
            self.dds_AOM_A1.sw.on()
        else:
            self.dds_AOM_A1.sw.off()
        delay(1 * ms)

        if self.AOM_A2_ON:
            self.dds_AOM_A2.sw.on()
        else:
            self.dds_AOM_A2.sw.off()
        delay(1 * ms)

        if self.AOM_A3_ON:
            self.dds_AOM_A3.sw.on()
        else:
            self.dds_AOM_A3.sw.off()
        delay(1 * ms)

        if self.AOM_A4_ON:
            self.dds_AOM_A4.sw.on()
        else:
            self.dds_AOM_A4.sw.off()
        delay(1 * ms)

        if self.AOM_A5_ON:
            self.dds_AOM_A5.sw.on()
        else:
            self.dds_AOM_A5.sw.off()
        delay(1 * ms)

        if self.AOM_A6_ON:
            self.dds_AOM_A6.sw.on()
        else:
            self.dds_AOM_A6.sw.off()
        delay(1 * ms)

        # D1: on Node2 this Urukul channel drives the GRIN2 AOM; the D1 light
        # itself is switched by ttl_D1_pumping_switch (gates, below). Unless
        # Node2_GRIN2_AOM_ON is ticked, the combined GRIN2 block at the end
        # turns this channel off again, so in the end the D1 box only opens
        # ttl_D1_pumping_switch, as in standalone bob.
        if self.D1_pumping_switch_ON:
            self.dds_GRIN2.sw.on()
        else:
            self.dds_GRIN2.sw.off()
        delay(1 * ms)

        # The combined GRIN1 block at the end decides how this channel ends.
        if self.GRIN1and2_ON:
            self.dds_GRIN1.sw.on()
        else:
            self.dds_GRIN1.sw.off()
        delay(1 * ms)

        ### TTL gates, active low: off() passes, on() blocks
        if self.Repump_AOM_switch_ON:
            self.ttl_repump_switch.off()
        else:
            self.ttl_repump_switch.on()
        delay(1 * ms)

        if self.pumping_repump_switch_ON:
            self.ttl_pumping_repump_switch.off()
        else:
            self.ttl_pumping_repump_switch.on()
        delay(1 * ms)

        if self.Excitation0_AOM_switch_ON:
            self.ttl_exc0_switch.off()
        else:
            self.ttl_exc0_switch.on()
        delay(1 * ms)

        # The combined GRIN1 block at the end decides how this gate ends.
        if self.GRIN1_AOM_switch_ON:
            self.ttl_GRIN1_switch.off()
        else:
            self.ttl_GRIN1_switch.on()
        delay(1 * ms)

        # The combined GRIN2 block at the end decides how this gate ends.
        if self.GRIN2_AOM_switch_ON:
            self.ttl_GRIN2_switch.off()
        else:
            self.ttl_GRIN2_switch.on()
        delay(1 * ms)

        if self.D1_pumping_switch_ON:
            self.ttl_D1_pumping_switch.off()  # D1 on
        else:
            self.ttl_D1_pumping_switch.on()  # D1 blocked
        delay(1 * ms)

        ### Combined GRIN modes (standalone bob). Each drives a GRIN DDS and
        ### its gate together. Ticked or not, it overrides the GRIN boxes
        ### above: unticked, that DDS ends off and that gate ends blocked.
        # dds_GRIN1 keeps the tone initialize_hardware() gave it,
        # f_GRIN1_D1_pumping_Node2 at p_GRIN1_D1_pumping_Node2.
        if self.Node2_GRIN1_AOM_ON:
            self.dds_GRIN1.sw.on()
            self.ttl_GRIN1_switch.off()
        else:
            self.dds_GRIN1.sw.off()
            self.ttl_GRIN1_switch.on()
        delay(1 * ms)

        if self.Node2_GRIN2_AOM_ON:
            # The one tone set here: GRIN2 runs at its excitation tone, not
            # the D1 tone initialize_hardware() gave this channel.
            self.dds_GRIN2.set(
                frequency=self.f_GRIN2_excitation_Node2,
                amplitude=dB_to_V(self.p_GRIN2_excitation_Node2),
            )
            self.dds_GRIN2.sw.on()
            self.ttl_GRIN2_switch.off()
        else:
            self.dds_GRIN2.sw.off()
            self.ttl_GRIN2_switch.on()
        delay(1 * ms)

    @kernel
    def set_microwaves(self):
        """Microwave and MW RF DDSes; both also need the safety box."""
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
    def set_coils(self):
        """MOT coils at Node2's persistent MOT calibration, or 0 V."""
        if self.disable_coils:
            self.zotino0.set_dac(
                [0.0, 0.0, 0.0, 0.0], channels=self.coil_channels
            )
        else:
            self.zotino0.set_dac(
                [
                    self.AZ_bottom_volts_MOT,
                    self.AZ_top_volts_MOT,
                    self.AX_volts_MOT,
                    self.AY_volts_MOT,
                ],
                channels=self.coil_channels,
            )
        delay(1 * ms)

    @kernel
    def apply_manual_state(self):
        """Initialize Node2 and set every output from the boxes."""
        # reset_core=False, as before. The run still resets the core:
        # run_feedback() and k10cr1_operations() start with core.reset().
        self.base.initialize_hardware(
            turn_off_dds_channels=True,
            turn_off_zotinos=True,
            initialize_spcms=False,
            reset_core=False,
        )
        self.everything_off()
        self.turn_on_AOMs()
        self.set_microwaves()
        self.set_coils()

    @kernel
    def run_feedback(self):
        """Feed back (or monitor) when the six fiber AOMs and the cooling DP
        are all requested on, then re-apply the boxes, combined GRIN modes
        included: the stabilizer switches AOMs while it measures. Re-applying
        sets none of the stabilized DDSes, so the powers the feedback
        corrected stay.
        """
        self.core.reset()
        if self._all_stabilized_aoms_requested:
            if self.run_laser_feedback:
                self.laser_stabilizer.run()
                delay(1 * ms)
                self.turn_on_AOMs()
            else:
                delay(10 * ms)
                self.laser_stabilizer.monitor()
                delay(1 * ms)
                self.turn_on_AOMs()

    @kernel
    def k10cr1_operations(self):
        """Rotate Node2's waveplates. Unlike the standalone utility, the
        852 target moves are gated on the 852 boxes (it used the 780 ones).
        """
        self.core.reset()
        delay(10 * ms)
        if self.go_to_home_780HWP:
            go_to_home(self, "780_HWP_Node2")
        if self.go_to_home_780QWP:
            go_to_home(self, "780_QWP_Node2")

        if self.go_to_target_780HWP:
            move_to_target_deg(
                self, name="780_HWP_Node2", target_deg=self.target_780_HWP
            )
        if self.go_to_target_780QWP:
            move_to_target_deg(
                self, name="780_QWP_Node2", target_deg=self.target_780_QWP
            )

        if self.move_780HWP_by:
            move_by_deg(
                self, name="780_HWP_Node2", target_deg=self.move_780_HWP_by
            )
        if self.move_780QWP_by:
            move_by_deg(
                self, name="780_QWP_Node2", target_deg=self.move_780_QWP_by
            )

        if self.go_to_home_852HWP:
            go_to_home(self, "852_HWP_Node2")
        if self.go_to_home_852QWP:
            go_to_home(self, "852_QWP_Node2")

        if self.go_to_target_852HWP:
            move_to_target_deg(
                self, name="852_HWP_Node2", target_deg=self.target_852_HWP
            )
        if self.go_to_target_852QWP:
            move_to_target_deg(
                self, name="852_QWP_Node2", target_deg=self.target_852_QWP
            )

        if self.move_852HWP_by:
            move_by_deg(
                self, name="852_HWP_Node2", target_deg=self.move_852_HWP_by
            )
        if self.move_852QWP_by:
            move_by_deg(
                self, name="852_QWP_Node2", target_deg=self.move_852_QWP_by
            )

        if self.go_to_optimized_852_settings:
            move_to_target_deg(
                self, name="852_HWP_Node2", target_deg=self.best_852HWP_to_max
            )
            move_to_target_deg(
                self, name="852_QWP_Node2", target_deg=self.best_852QWP_to_max
            )

    def run(self):
        self.apply_manual_state()
        self.run_feedback()
        if self._k10cr1_requested:
            self.k10cr1_operations()
