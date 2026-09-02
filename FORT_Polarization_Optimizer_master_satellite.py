"""FORT polarization optimization for the master-satellite architecture.

Physics (unchanged from the standalone FORT_Polarization_Optimizer):
a few percent of the FORT light passes the parabolic mirror through a small
hole, then a polarizer pre-aligned to the preferred FORT polarization, and is
collected in the MM fiber. Finding the 852 HWP/QWP angles that maximize the
MM power optimizes the polarization. The scan is a shrinking zigzag grid; the
MM power is normalized by the FORT APD power so the metric depends on
polarization, not overall power drift.

which_node selects Node1 or Node2. The selected node's 852 waveplates are
scanned through the single k10cr1_ndsp controller (node-suffixed axis names),
powers are read from the node's own samplers via the unchanged legacy helpers
(the dataset-redirect layer stores every monitor and calibration under the
node-suffixed names: FORT_MM_monitor_NodeX, best_852HWP_to_max_NodeX, ...).

FUTURE REFERENCE - two_nodes design (deliberately not implemented yet;
revisit after two-node PARALLEL FORT laser feedback exists, see
plan_codex_detail.md section 18):

- which_node gains "two_nodes" and Base is configured in two_nodes mode.
- Both nodes march through one grid by point index. Per scan point, all four
  waveplates (852_HWP_Node1, 852_QWP_Node1, 852_HWP_Node2, 852_QWP_Node2)
  are set FIRST, then each node records its own powers. The k10cr1 moves are
  blocking RPCs on the single controller, so sequential moves already
  guarantee all-set-before-record; truly concurrent rotation needs an async
  move method in ndsp/k10cr1/k10cr1_driver.py plus a server restart.
- Per-node measurement devices (suffixed presentation):
  MM  = sampler1_NodeX channel 7 (both nodes);
  APD = sampler1_Node1 / sampler0_Node2, channel 6 (the legacy alice/bob
  asymmetry). The legacy record helpers hard-code self.sampler0/1 and
  which_node, so two_nodes mode needs parameterized kernels such as
  _record_power(sampler, channel) instead.
- Each node normalizes MM by (APD / set_point_FORT_APD_loading_NodeX) and
  tracks its own best; after the first iteration each node's grid re-centers
  on its own best (same shape and step sizes, same point index).
- The per-row early-stop shortcut (count_down/stop_loop) is per-node state;
  either disable it in two_nodes or stop a row only when every node's
  criterion fires.
- Run both nodes' FORT feedback in parallel before/during the scan once
  section 18 allows it; until then two_nodes would have to run at the
  persistent default FORT powers.
- Finish per node: move its plates to its own best and persist
  best_852HWP_to_max_NodeX, best_852QWP_to_max_NodeX, best_852_power_NodeX,
  best_852_power_ref_NodeX.
"""

from artiq.experiment import *

import numpy as np

from subroutines.k10cr1_functions import (
    go_to_home,
    move_to_target_deg,
    record_FORT_APD_power,
    record_FORT_MM_power,
    time_to_rotate_in_ms,
)
from utilities.BaseExperiment_master_satellite import (
    BaseExperimentMasterSatellite,
    _DatasetRedirectMixin,
)


class FORT_Polarization_Optimizer_master_satellite(
    _DatasetRedirectMixin, EnvExperiment
):
    """FORT_Polarization_Optimizer_master_satellite

    Optimize the selected node's FORT polarization via its 852 waveplates.
    """

    NODE_BY_SELECTION = {"node1": "Node1", "node2": "Node2"}
    # Tests inject a fake stabilizer class here; None selects the real
    # AOMPowerStabilizer from subroutines/aom_feedback.py.
    _stabilizer_factory = None

    def build(self):
        # Repository examination supplies None for arguments and may run with
        # an empty dataset namespace; bind the suffixed device superset now
        # and defer mode validation and dataset loading to prepare().
        self.base = BaseExperimentMasterSatellite(experiment=self)
        self.base.build()

        self.setattr_argument(
            "which_node",
            EnumerationValue(tuple(self.NODE_BY_SELECTION)),
            "Node selection",
        )
        self.setattr_argument(
            "initialize_to_home", BooleanValue(False), "Initialization"
        )
        self.setattr_argument(
            "tolerance_deg",
            NumberValue(1, ndecimals=2, step=1),
            "optimization parameters",
        )
        self.setattr_argument(
            "full_range",
            NumberValue(90, ndecimals=1, step=1),
            "optimization parameters",
        )
        self.setattr_argument(
            "sample_pts",
            NumberValue(7, ndecimals=1, step=1),
            "optimization parameters",
        )
        self.setattr_argument(
            "search_start_from_scratch",
            BooleanValue(False),
            "optimization settings",
        )

    def prepare(self):
        # tolerance_deg/full_range/sample_pts are run-local GUI values that
        # share their names with persistent node datasets. Capture the
        # submitted values before configure_execution loads those datasets
        # over the same attributes, and re-assert them afterward so the GUI
        # wins for this run.
        submitted_tolerance_deg = self.tolerance_deg
        submitted_full_range = self.full_range
        submitted_sample_pts = self.sample_pts

        try:
            node = self.NODE_BY_SELECTION[str(self.which_node)]
        except KeyError:
            raise ValueError(
                f"Unsupported which_node {self.which_node!r}; expected "
                "'node1' or 'node2'. two_nodes is not implemented yet - see "
                "the future-reference notes in this file."
            ) from None

        self.base.configure_execution("single_node", node)

        # The legacy record helpers branch on the alice/bob presentation.
        self.which_node = self.base.NODE_LEGACY_NAMES[node]
        self.base.prepare()
        self.base.prepare_laser_stabilizer(
            stabilizer_factory=self._stabilizer_factory
        )

        # Re-assert after base.prepare(): these are node variables, so the
        # compatibility refresh inside prepare() re-projects the dataset
        # values over them.
        self.tolerance_deg = submitted_tolerance_deg
        self.full_range = submitted_full_range
        self.sample_pts = submitted_sample_pts

        # Base publishes k10cr1_ndsp as a lazy node-aware proxy (it connects
        # only on first use). Axis names stay explicitly suffixed here; the
        # proxy passes already-suffixed names through unchanged.
        self._axis_852_HWP = f"852_HWP_{node}"
        self._axis_852_QWP = f"852_QWP_{node}"

    @kernel
    def initialize_hardware(self):
        self.base.initialize_hardware()

        if self.initialize_to_home:
            go_to_home(self, self._axis_852_HWP)
            go_to_home(self, self._axis_852_QWP)
            print("Both Waveplates Initialized to Home")
        else:
            #### initializing the waveplates to best_852 parameters is
            #### necessary if using MM feedback; this will do nothing if it
            #### has been optimized before
            move_to_target_deg(
                self, name=self._axis_852_HWP,
                target_deg=self.best_852HWP_to_max,
            )
            move_to_target_deg(
                self, name=self._axis_852_QWP,
                target_deg=self.best_852QWP_to_max,
            )

    def initialize_datasets(self):
        # Written through the dataset-redirect layer, so these land under the
        # node-suffixed names (FORT_MM_monitor_NodeX, ...).
        self.set_dataset("FORT_MM_monitor", [], broadcast=True)
        self.set_dataset("FORT_APD_monitor", [], broadcast=True)
        self.set_dataset("HWP_angle", [], broadcast=True)
        self.set_dataset("QWP_angle", [], broadcast=True)

    @kernel
    def optimization_routine_zigzag_power_normalized(self):
        """Shrinking zigzag grid scan on normalized FORT MM power.

        Ported from the standalone optimizer: as the routine proceeds to the
        next HWP value it reverses the QWP scan direction, minimizing total
        rotation time. The MM power is normalized by
        power_APD/set_point_FORT_APD_loading so the metric depends mainly on
        polarization, not overall power drift.
        """

        self.core.reset()
        delay(1 * s)

        # run stabilizer
        self.core.break_realtime()
        if self.enable_laser_feedback:
            self.stabilizer_FORT.run(setpoint_index=0)  # FORT loading setpoint

        ### turning FORT on in the beginning to give enough time to stabilize.
        self.dds_FORT.set(
            frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude
        )
        delay(100 * us)
        self.dds_FORT.sw.on()  ### turns FORT on

        delay(1 * ms)
        self.dds_cooling_DP.sw.off()  ### turn off cooling
        self.ttl_repump_switch.on()  ### turn off MOT RP

        delay(1 * ms)
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()
        delay(1 * ms)

        tolerance = float(self.tolerance_deg)  # stop when step < tolerance
        full_range = float(self.full_range)  # Start with a full range
        sample_pts = int(self.sample_pts)  # Number of samples per iteration

        # only half_range used in iteration. do not use full_range
        half_range = full_range / 2  # Half the full range
        measurement_count = 0  # Track number of power measurements

        if self.search_start_from_scratch:
            best_HWP, best_QWP, best_power = 0.0, 0.0, 0.0
        else:
            best_HWP, best_QWP, best_power = (
                self.best_852HWP_to_max, self.best_852QWP_to_max, 0.0
            )
        previous_hwp, previous_qwp, previous_power = best_HWP, best_QWP, 0.0

        power = 0.0
        power_APD = 0.0

        tolerance_satisfied = False
        iteration = 0

        # Iterative search loop
        while not tolerance_satisfied:
            print("starting iteration no.", iteration)
            delay(10 * ms)
            time_ite = time_to_rotate_in_ms(self, half_range)

            delay(2 * time_ite * ms + 1 * s)  # rotate two waveplates

            steps = half_range * 2 / (sample_pts - 1)

            # power difference with steps less than 1 is negligible.
            if steps < self.tolerance_deg:
                tolerance_satisfied = True

                steps = self.tolerance_deg
                # how many sample pts to cover the range
                sample_pts = int(half_range / steps) * 2 + 1
                half_range = steps * (sample_pts - 1) / 2

            hwp_values = np.array(
                [best_HWP - half_range + i * steps for i in range(sample_pts)]
            )
            qwp_values = np.array(
                [best_QWP - half_range + i * steps for i in range(sample_pts)]
            )

            for hwp_i in range(len(hwp_values)):
                # time to rotate hwp
                current_hwp = hwp_values[hwp_i]
                time_hwp = time_to_rotate_in_ms(self, current_hwp - previous_hwp)
                delay(time_hwp * ms)
                count_down = 0

                stop_loop = False

                for qwp_i in range(len(qwp_values)):

                    if hwp_i % 2 != 0:
                        qwp_i = len(qwp_values) - 1 - qwp_i
                    if not stop_loop:
                        current_qwp = qwp_values[qwp_i]
                        # time to rotate qwp
                        time_qwp = time_to_rotate_in_ms(
                            self, current_qwp - previous_qwp
                        )
                        delay(time_qwp * ms)

                        self.append_to_dataset("HWP_angle", current_hwp)
                        self.append_to_dataset("QWP_angle", current_qwp)

                        with parallel:
                            move_to_target_deg(
                                self, name=self._axis_852_HWP,
                                target_deg=current_hwp,
                            )
                            move_to_target_deg(
                                self, name=self._axis_852_QWP,
                                target_deg=current_qwp,
                            )

                        delay(1 * s)

                        power = record_FORT_MM_power(self)
                        power_APD = record_FORT_APD_power(self)

                        power_APD_norm = (
                            power_APD / self.set_point_FORT_APD_loading
                        )

                        power = power / power_APD_norm

                        if half_range > 5:
                            delay(1 * s)

                        if power > best_power:  # Update best if improved
                            best_power = power
                            best_HWP = current_hwp
                            best_QWP = current_qwp

                        previous_qwp = current_qwp
                        measurement_count += 1  # Increment measurement counter

                        if power < previous_power:
                            count_down += 1
                            if count_down > int(sample_pts / 2):
                                stop_loop = True

                        previous_power = power

                previous_hwp = current_hwp
                delay(1 * s)

            half_range = steps

            delay(1 * s)
            print(
                "iteration # ", iteration,
                " : best_HWP, best_QWP, best_power = ",
                best_HWP, ", ", best_QWP, ", ", best_power,
            )
            delay(10 * ms)
            iteration += 1

        # move back to the best HWP, QWP
        move_to_target_deg(
            self, name=self._axis_852_HWP, target_deg=best_HWP
        )
        move_to_target_deg(
            self, name=self._axis_852_QWP, target_deg=best_QWP
        )

        print(
            "previous best_HWP, best_QWP, best_power = ",
            self.best_852HWP_to_max, ", ", self.best_852QWP_to_max, ", ",
            self.best_852_power,
        )
        delay(10 * ms)

        self.core.break_realtime()
        self.dds_FORT.sw.off()  ### turns FORT off
        delay(1 * ms)

        # Persisted through the dataset-redirect layer to the node-suffixed
        # calibration datasets (best_852HWP_to_max_NodeX, ...).
        self.set_dataset(
            "best_852HWP_to_max", best_HWP, broadcast=True, persist=True
        )
        self.set_dataset(
            "best_852QWP_to_max", best_QWP, broadcast=True, persist=True
        )
        self.set_dataset(
            "best_852_power", best_power, broadcast=True, persist=True
        )

    @kernel
    def run_feedback_and_record_ref_power(self):
        """Stabilize the FORT to the loading setpoint, then record the
        reference powers at the freshly optimized waveplate angles."""
        delay(1 * s)
        self.core.reset()

        delay(0.1 * ms)
        self.stabilizer_FORT.run(setpoint_index=0)  # FORT loading setpoint
        delay(0.1 * ms)
        self.dds_FORT.set(
            frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude
        )
        self.dds_FORT.sw.on()  ### turns FORT on
        delay(0.1 * ms)

        power = record_FORT_MM_power(self)
        record_FORT_APD_power(self)

        self.dds_FORT.sw.off()  ### turns FORT off
        delay(0.1 * ms)

        print("After feedback - best_852_power set to ", power)
        delay(10 * ms)
        self.set_dataset(
            "best_852_power_ref", power, broadcast=True, persist=True
        )

    def run(self):
        # Unlike the standalone optimizer, prepare() is not re-run here: the
        # master-satellite Base prepare is one-shot, and initialize_hardware
        # performs the core reset.
        self.initialize_hardware()
        self.initialize_datasets()

        self.optimization_routine_zigzag_power_normalized()

        self.run_feedback_and_record_ref_power()  # recording reference power
