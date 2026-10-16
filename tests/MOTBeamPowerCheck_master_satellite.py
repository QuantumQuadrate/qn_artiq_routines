"""
Sample the six MOT beam photodiodes (AOM A1-A6) of one node, one beam at a
time, at the powers the last feedback run left -- without changing anything.

What it does, per round:
  * cooling double-pass on, repump RF blocked, every MOT AOM off;
  * for A1, A2, ... A6 in turn: that AOM on, wait its t_measure_delay, sample,
    record the raw photodiode voltage, AOM off.
This is exactly how AOMPowerStabilizer measures its "series" channels
(subroutines/aom_feedback.py, AOMPowerStabilizer.run), minus the feedback
step: every AOM runs at its stored frequency and p_AOM_A*_<node> power, which
is where feedback leaves it, and nothing is written back. The powers are
re-read when the run STARTS, not when it is prepared: ARTIQ prepares a queued
run while the previous one is still running, so a prepare-time value could
predate that run's feedback. The dB values applied are recorded in
MOT_beam_p_AOM_dB_<node> and printed with the readings.

Which Sampler card and channel each beam is read on comes from the node's
feedback config (utilities/config/<alice|bob>/feedback_channels.json) through
the same stabilizer the experiments use, so a rewiring is picked up
automatically. prepare() prints each beam's card and channel.

Results: MOT_beam_PD_A1_<node> ... MOT_beam_PD_A6_<node>, raw volts, one
value per round, all six in ONE applet ("MOT beam PDs <node>", group "Both
nodes", created only when show_applet is ticked). analyze() prints each
beam's mean, spread and set point.

At the end the node is left in the MOT-light state feedback leaves it in:
cooling DP and the six MOT AOMs on, repump RF enabled. The other node is not
touched.
"""

from artiq.experiment import *
import logging
import numpy as np

import sys, os
### ARTIQ puts only this file's own directory on sys.path; the worker's cwd is
### artiq-master, so add the repository from there before any repo import.
cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd + "\\repository\\qn_artiq_routines")

from utilities.BaseExperiment_master_satellite import (
    BaseExperimentMasterSatellite,
    _DatasetRedirectMixin,
)
from utilities.conversions import dB_to_V


class MOTBeamPowerCheck_master_satellite(_DatasetRedirectMixin, EnvExperiment):
    """MOTBeamPowerCheck_master_satellite
    Sample one node's six MOT beam photodiodes one beam at a time, at the
    powers feedback left, without changing them.
    """

    VALID_NODES = ("Node1", "Node2")
    BEAMS = ("dds_AOM_A1", "dds_AOM_A2", "dds_AOM_A3",
             "dds_AOM_A4", "dds_AOM_A5", "dds_AOM_A6")

    def build(self):
        self.base = BaseExperimentMasterSatellite(experiment=self)
        self.base.build()

        self.setattr_argument(
            "selected_node",
            EnumerationValue(self.VALID_NODES, default="Node2"),
            "Node selection",
            tooltip="Node1 = alice, Node2 = bob. Only this node's hardware is "
                    "initialized and used; the other node is left alone.",
        )
        self.setattr_argument(
            "n_rounds",
            NumberValue(100, type='int', scale=1, ndecimals=0, step=1),
            tooltip="How many times to sample all six beams.",
        )
        self.setattr_argument(
            "t_between_rounds_ms",
            NumberValue(10, type='int', scale=1, ndecimals=0, step=1),
            tooltip="Pause after each round of six readings.",
        )
        self.setattr_argument(
            "show_applet",
            BooleanValue(True),
            "Applets",
            tooltip="Open the 'MOT beam PDs <node>' applet with all six "
                    "beams. Requires the applet dock's CCB policy to be "
                    "'Create and enable/disable applets'.",
        )
        print("build - done")

    def prepare(self):
        node = str(self.selected_node)
        if node not in self.VALID_NODES:
            raise ValueError(
                f"Unsupported selected_node {self.selected_node!r}; expected "
                "'Node1' or 'Node2'."
            )
        self.base.configure_execution("single_node", node)
        ### The shared feedback code reads the alice/bob presentation.
        self.which_node = self.base.NODE_LEGACY_NAMES[node]
        self.base.prepare()
        ### Builds self.laser_stabilizer from the node's feedback config, with
        ### every channel at its stored frequency and p_AOM_A*_<node> power.
        self.base.prepare_laser_stabilizer()

        channels = {ch.name: ch for ch in self.laser_stabilizer.all_channels}
        missing = [name for name in self.BEAMS if name not in channels]
        if missing:
            raise RuntimeError(
                f"{node}'s laser stabilizer has no channel for {missing}. "
                f"They come from fast_feedback_dds_list_{node} and "
                f"utilities/config/{self.which_node}/feedback_channels.json; "
                "both must list every MOT AOM."
            )
        self.beams = [channels[name] for name in self.BEAMS]

        self._node = node
        self.beam_datasets = [
            f"MOT_beam_PD_A{k}_{node}" for k in range(1, len(self.BEAMS) + 1)
        ]
        self.labels_dataset = "MOT_beam_PD_labels"
        self.applied_power_dataset = f"MOT_beam_p_AOM_dB_{node}"
        self.applied_dB = []
        self.points_dataset = "MOT_beam_PD_points"
        self.beam_set_points = [float(ch.set_points[0]) for ch in self.beams]

        ### Seeded on the host, never in the kernel: append_to_dataset can only
        ### mutate a dataset that already exists in this worker, and
        ### set_dataset from a kernel is a blocking RPC that spends slack. The
        ### applet drops each series' first element, the [0.0] seed.
        for name in self.beam_datasets:
            self.set_dataset(name, [0.0], broadcast=True)
        self.set_dataset(self.labels_dataset,
                         [f"AOM{k}" for k in range(1, len(self.BEAMS) + 1)],
                         broadcast=True)
        self.set_dataset(self.points_dataset, int(self.n_rounds), broadcast=True)

        if self.show_applet:
            self._create_applet()

        for ch, set_point in zip(self.beams, self.beam_set_points):
            sampler = ch.buffer_index // 8
            print(f"{node} {ch.name}: stabilizer sampler #{sampler}, "
                  f"ch {ch.buffer_index % 8}, set point {set_point}")
        print("prepare - done")

    def _create_applet(self):
        """One applet, all six beams. Never fatal: it is a convenience, and
        artiq_run / the offline compile check have no real dashboard."""
        try:
            from applets_master_satellite import (
                AppletSpec,
                create_shared_applets_for,
            )
            spec = AppletSpec(
                f"MOT beam PDs {self._node}",
                "plot_xy_multichannel.py",
                (self.beam_datasets[0], self.points_dataset),
                tuple((f"--y{k}", name) for k, name in
                      enumerate(self.beam_datasets[1:], start=2))
                + (("--labels", self.labels_dataset),),
                None,
            )
            ### The shared group is never retired by other experiments, so
            ### the applet stays where it is until closed by hand.
            create_shared_applets_for(self, self.base, shared_specs=(spec,))
        except Exception as error:  # noqa: BLE001 - convenience only
            logging.warning("could not create the applet: %s", error)

    def run(self):
        ### The powers feedback left, read NOW rather than at prepare: this run
        ### may have been prepared while an earlier one was still running
        ### feedback. ch.dB_dataset is the name feedback itself writes
        ### (AOMPowerStabilizer.update_dB_dataset), and the dataset redirect
        ### resolves it to this node's p_AOM_A*_<node> both ways.
        self.applied_dB = []
        for ch in self.beams:
            dB = float(self.get_dataset(ch.dB_dataset, archive=False))
            ch.amplitude = dB_to_V(dB)
            ch.amplitudes[0] = ch.amplitude
            self.applied_dB.append(dB)
        ### A name nothing else writes, so recording it locally is safe.
        self.set_dataset(self.applied_power_dataset, self.applied_dB)
        self.measure_beams()

    @kernel
    def measure_beams(self):
        ### Owns the core reset; waits for the satellite when Node2 is
        ### selected; sets every DDS to its stored frequency and power, off.
        self.base.initialize_hardware()
        self.core.break_realtime()
        delay(1 * ms)

        ### As AOMPowerStabilizer.run does before a series measurement: the
        ### photodiodes must see the cooling light only.
        self.ttl_repump_switch.on()            # block RF to the repump AOM
        self.ttl_pumping_repump_switch.on()    # block RF to the pumping repump AOM
        delay(1 * ms)
        for ch in self.beams:
            ch.set_dds_to_defaults()
            delay(0.1 * ms)
            ch.dds_obj.sw.off()
            delay(0.1 * ms)
        self.dds_cooling_DP.sw.on()
        delay(1 * ms)

        for n in range(self.n_rounds):
            ### Each round spends real time on six samplings and six async
            ### appends; re-arm before the first output of the next round.
            self.core.break_realtime()
            k = 0
            for ch in self.beams:
                ch.dds_obj.sw.on()
                delay(ch.t_measure_delay)      # detector rise time
                self.laser_stabilizer.measure()
                delay(0.1 * ms)
                self.append_to_dataset(
                    self.beam_datasets[k],
                    self.laser_stabilizer.measurement_array[ch.buffer_index])
                ch.dds_obj.sw.off()
                delay(0.1 * ms)
                k += 1
            delay(self.t_between_rounds_ms * ms)

        ### Leave the MOT light on, as feedback does (leave_MOT_AOMs_on).
        self.core.break_realtime()
        self.ttl_repump_switch.off()           # enable RF to the repump AOM
        delay(0.1 * ms)
        for ch in self.beams:
            ch.dds_obj.sw.on()
            delay(10 * us)
        self.dds_cooling_DP.sw.on()
        delay(1 * ms)

    def analyze(self):
        print(f"MOT beam photodiodes, {self._node}, raw volts:")
        applied = self.applied_dB or [None] * len(self.beams)
        for name, ch, set_point, dB in zip(self.beam_datasets, self.beams,
                                           self.beam_set_points, applied):
            power = "" if dB is None else f"p {dB:.3f} dB, "
            values = np.asarray(self.get_dataset(name)[1:], dtype=float)
            if values.size == 0:
                print(f"  {ch.name}: {power}no readings")
                continue
            print(f"  {ch.name}: {power}mean {values.mean():.4f} V, "
                  f"std {values.std():.4f} V, n {values.size}, "
                  f"set point {set_point:.4f} V, "
                  f"mean/set point {values.mean() / set_point:.3f}")
