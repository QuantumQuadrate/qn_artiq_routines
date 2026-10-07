"""Two-node experiments for the master-satellite system.

HOW TO READ THIS FILE
---------------------
Same shape as section 4 of experiment_functions.py, deliberately. A two-node
run has always been "the same code running on both nodes"; under one DRTIO
kernel that becomes two explicit lines, and EVERY DEVICE NAMES ITS NODE:

    self.dds_AOM_A1_Node1.sw.off()      turn A1 off on Node1
    self.dds_AOM_A1_Node2.sw.off()      ... and on Node2

There is no bare `self.dds_AOM_A1` in two-node mode. That is the whole
convention, and it is enforced rather than remembered: Base publishes no bare
device aliases here, so a bare device name fails to COMPILE instead of
resolving to something.

An earlier version of this file did have bare names, published as broadcast
objects that fanned out to both crates. It read more like the legacy code, and
it cost more than it saved: a bare name became LEGAL everywhere in the stack,
including in shared modules. That is exactly how laser_stabilizer_Node2.run()
came to drive NODE1's repump -- aom_feedback.py reached for the bare
self.exp.ttl_repump_switch and got a broadcast object instead of an error, and
the bug surfaced only because an unrelated RTIOUnderflow happened to print the
stack. Two sequential writes land at the same now_mu anyway, so the explicit
pair is not just safer, it is the same hardware behaviour.

The four SPCMs are the one deliberate exception, and they are not a broadcast:
all four are master-local and all four already see both nodes' fluorescence
through the beamsplitter fan-out, so ttl_SPCM0_counter is the single real
counter, exactly as in single-node mode.

THREE PLACES TWO NODES CANNOT SHARE ONE LINE AT ALL
---------------------------------------------------
Beyond naming, these need genuinely different code per node:

1. VALUES. The nodes disagree about nearly every setting that matters --
   f_FORT 245 vs 240 MHz, f_cooling_DP_RO 120.37 vs 130 MHz, every coil
   voltage, every ampl_* (they derive from per-node dBm calibrations).

2. DURATIONS. t_PGC_after_loading is 0.6 ms on Node1 and 1.0 ms on Node2;
   t_recooling_after_first_shot is 0.0 and 1.0 ms; t_FORT_drop is per node
   too. The two nodes' sequences are different LENGTHS. Two independent
   kernels could each just delay() its own amount; one kernel cannot. Those
   stages are placed from a common t0 with at_mu and then resynced, which is
   the only place the code stops being a single linear list of steps.

3. BRANCHES on a per-node flag. do_PGC_after_loading gates a stage whose
   duration already differs, and PGC_and_RO_with_on_chip_beams selects A5/A6,
   which are a per-node AOM pair. Each node's branch is written out, so one
   node can take a stage the other skips.

SCALARS the two nodes agree on are still published bare by Base, and that is
about values rather than devices: t_SPCM_first_shot and friends are consumed by
ONE hardware event (a single gate of all four master-local SPCMs, or a single
interval on the one shared timeline), so there is only one value to have. Base
RAISES if the nodes ever diverge on one, rather than handing the sequence one
node's value, and the check runs after overrides and at every scan point --
not only at prepare time. A duration that merely happens to be equal today
does not qualify. See BaseExperimentMasterSatellite.TWO_NODE_AGREEING_SCALARS.

WHAT IS COPIED FROM experiment_functions.py
-------------------------------------------
first_shot, second_shot, the loader and end_measurement are all COPIED into
this file, not imported. The single-node originals cannot be called here: they
reach for self.coil_channels, bare self.ampl_*, and which_node-flavoured
branches, none of which exist in two-node mode. end_measurement is a near
copy -- only the dead magnetometer branch is dropped -- but it is still a copy,
so a fix to the original does NOT propagate here. Check both when changing
either.

FEEDBACK
--------
Works as it does in standalone, per node, writing the SAME per-node datasets
(p_AOM_A1_Node1, FORT_monitor_Node2, ...), so the applets and the analysis
notebooks are unchanged. Base builds one AOMPowerStabilizer per crate and
publishes laser_stabilizer_Node1/_Node2 and stabilizer_FORT_Node1/_Node2.
"""

from artiq.experiment import *

import numpy as np


MASTER_SATELLITE_SANITY_ATTRIBUTES = (
    "core",
    "dds_FORT_Node1",
    "dds_FORT_Node2",
    "sampler0_Node1",
    "sampler0_Node2",
    "zotino0_Node1",
    "zotino0_Node2",
    "SPCM_H1",
    "SPCM_V1",
    "SPCM_H2",
    "SPCM_V2",
    "SPCM_H1_counter",
    "SPCM_V1_counter",
    "SPCM_H2_counter",
    "SPCM_V2_counter",
    "f_FORT_Node1",
    "f_FORT_Node2",
    "ttl_repump_switch_Node1",
    "ttl_repump_switch_Node2",
)


def master_satellite_namespace_sanity_experiment(self):
    """Validate the initial two-node namespace without touching hardware.

    This deliberately performs only host-side attribute inspection. Device
    initialization, core reset, DRTIO readiness, TTL direction, and output
    state remain responsibilities of BaseExperimentMasterSatellite.
    """
    missing = [
        name
        for name in MASTER_SATELLITE_SANITY_ATTRIBUTES
        if not hasattr(self, name)
    ]
    if missing:
        raise RuntimeError(
            "Master-satellite namespace sanity check failed; missing "
            "attribute(s): " + ", ".join(missing)
        )
    return True


#####################################################

@kernel
def record_FORT_powers_both_nodes(self):
    """
    Record each node's FORT MM and APD powers into its own dataset.

    Two-node version of record_FORT_MM_power / record_FORT_APD_power from
    k10cr1_functions.py. Those cannot be imported: each reads one node's
    sampler and appends to ONE FORT_MM_monitor / FORT_APD_monitor, and the
    APD one branches on which_node to choose sampler1 vs sampler0. Here both
    nodes are measured into the per-node datasets a single-node run writes,
    so the FORT monitor applets keep working.

    Wiring, from the single-node functions: MM is sampler1 ch7 on BOTH nodes
    ("Node1 & Node2 both have FORT_MM connected at Sampler1. ch7"); APD is
    sampler1 ch6 on Node1 and sampler0 ch6 on Node2.
    """
    measurement_buf = np.array([0.0] * 8)
    avgs = 200

    ### FORT MM -- sampler1 ch7 on both nodes
    mm_node1 = 0.0
    mm_node2 = 0.0
    for i in range(avgs):
        self.sampler1_Node1.sample(measurement_buf)
        mm_node1 += measurement_buf[self.FORT_MM_sampler_ch_Node1]
        delay(0.1 * ms)
        self.sampler1_Node2.sample(measurement_buf)
        mm_node2 += measurement_buf[self.FORT_MM_sampler_ch_Node2]
        delay(0.1 * ms)
    self.append_to_dataset("FORT_MM_monitor_Node1", mm_node1 / avgs)
    self.append_to_dataset("FORT_MM_monitor_Node2", mm_node2 / avgs)
    delay(0.1 * ms)

    ### FORT APD -- ch6, on sampler1 for Node1 and sampler0 for Node2
    apd_node1 = 0.0
    apd_node2 = 0.0
    for i in range(avgs):
        self.sampler1_Node1.sample(measurement_buf)
        apd_node1 += measurement_buf[6]
        delay(0.1 * ms)
        self.sampler0_Node2.sample(measurement_buf)
        apd_node2 += measurement_buf[6]
        delay(0.1 * ms)
    self.append_to_dataset("FORT_APD_monitor_Node1", apd_node1 / avgs)
    self.append_to_dataset("FORT_APD_monitor_Node2", apd_node2 / avgs)
    delay(0.1 * ms)

    self.core.break_realtime()


@kernel
def run_feedback_and_record_FORT_MM_power(self, record_power=True):
    """
    Function:
        1. Runs feedback to everything in list (6 MOT AOMs and 3 FORT setpoints), on both nodes
        2. records FORT MM and APD powers

    * IF you want to only run feedback and disable recording, set "record_power = False"
    """
    ### EVERY run() gets its own break_realtime, not a fixed delay.
    ###
    ### AOMPowerStabilizer.run() is wall-clock-heavy host work -- get_dataset,
    ### set_dataset, append_to_dataset and print_async are all host round trips
    ### -- while its coded timeline advances only by the delays written inside
    ### it. Wall clock runs ahead of the cursor, so slack is spent and never
    ### recovered. Six runs in a row, separated by nothing but 0.1 ms delays,
    ### is how the 2026-10-07 two-node loading run hit -5.698 ms on channel 5:
    ### the first RTIO event of laser_stabilizer_Node2.run() -- the repump
    ### switch at aom_feedback.py:550 -- was submitted ~5.7 ms in the past.
    ### Node1's run and Node2's run had no delay between them at all.
    ###
    ### A delay() cannot fix this: it would have to be larger than the host
    ### work it is guessing at, and the overrun scales with dataset size and
    ### master load. break_realtime resyncs the cursor to the wall clock, so
    ### each run starts with full slack regardless of what the previous one
    ### cost. Feedback is not timing-critical relative to the sequence -- it
    ### runs before the measurement loop, or between loading attempts -- so
    ### breaking timeline continuity here costs nothing.
    ### BOTH nodes at BOTH setpoints. FeedbackChannel.amplitudes is
    ### np.zeros(len(set_points)) with only [0] seeded (aom_feedback.py:106),
    ### and run(setpoint_index=N) writes ONLY amplitudes[N]. So an index that
    ### is never fed back stays 0.0 -- which, read back as a dds amplitude,
    ### turns the FORT OFF rather than failing.
    ###
    ### The single-node original (experiment_functions.py:413-417) runs the one
    ### FORT stabilizer twice, at index 2 then index 1, then laser_stabilizer
    ### at index 0, so all three are populated. Interleaving two nodes must
    ### keep that second axis: run each node at 2 AND at 1. Pairing one node
    ### with one index instead -- Node1 always 2, Node2 always 1 -- leaves
    ### stabilizer_FORT_Node1.amplitudes[1] at 0.0, and first_shot,
    ### second_shot and two_node_alternating_shot all read exactly that index,
    ### so Node1's FORT would be commanded to zero amplitude at every readout
    ### and drop its atom. Loading is unaffected: it reads .amplitude, which
    ### is amplitudes[0].
    ###
    ### laser_stabilizer runs stay LAST so the FORT is left at the loading
    ### setpoint, as in the single-node original.
    self.core.break_realtime()
    self.stabilizer_FORT_Node1.run(setpoint_index=2)  # FORT science holding setpoint
    self.core.break_realtime()
    self.stabilizer_FORT_Node2.run(setpoint_index=2)  # FORT science holding setpoint
    self.core.break_realtime()
    self.stabilizer_FORT_Node1.run(setpoint_index=1)  # FORT science setpoint
    self.core.break_realtime()
    self.stabilizer_FORT_Node2.run(setpoint_index=1)  # FORT science setpoint
    self.core.break_realtime()
    self.laser_stabilizer_Node1.run()  # 6 MOT AOMs and FORT loading setpoint
    self.core.break_realtime()
    self.laser_stabilizer_Node2.run()  # 6 MOT AOMs and FORT loading setpoint
    self.core.break_realtime()

    ## if laser_stabilizer.run() is in the last sequence, it will leave the FORT at loading setpoint.

    ## record FORT MM and APD powers
    if record_power:
        self.dds_FORT_Node1.sw.on()  ### turns FORT on, both nodes
        self.dds_FORT_Node2.sw.on()
        delay(0.1*ms)

        record_FORT_powers_both_nodes(self)

    ### record_FORT_powers_both_nodes does its own host RPCs, and with
    ### record_power=False the last resync above is already several RPCs old,
    ### so restore slack for the CALLER's next RTIO event. Covers every call
    ### site at once.
    self.core.break_realtime()


@kernel
def first_shot(self):
    """
    first atom readout, both nodes at once.

    Copied from experiment_functions.first_shot. The switching lines are
    unchanged -- the bare names drive both crates. Only the two .set() calls
    are doubled, because f_FORT and f_cooling_DP_RO differ per node.

    Turns on:  Cooling DP, MOT RP, all 6 fiber AOMs (both nodes)
    Turns off at the end:  Cooling DP, MOT RP
    """

    ### set the FORT AOM to the science settings
    self.dds_FORT_Node1.set(frequency=self.f_FORT_Node1,
                            amplitude=self.stabilizer_FORT_Node1.amplitudes[1])
    self.dds_FORT_Node2.set(frequency=self.f_FORT_Node2,
                            amplitude=self.stabilizer_FORT_Node2.amplitudes[1])

    ### set the cooling DP AOM to the readout settings
    self.dds_cooling_DP_Node1.set(frequency=self.f_cooling_DP_RO_Node1,
                                  amplitude=self.ampl_cooling_DP_RO_Node1)
    self.dds_cooling_DP_Node2.set(frequency=self.f_cooling_DP_RO_Node2,
                                  amplitude=self.ampl_cooling_DP_RO_Node2)
    delay(5 * us)

    self.ttl_repump_switch_Node1.off()  ### turn on MOT RP
    self.ttl_repump_switch_Node2.off()
    self.dds_cooling_DP_Node1.sw.on()  ### Turn on cooling
    self.dds_cooling_DP_Node2.sw.on()
    delay(5 * us)
    delay(0.1 * ms)

    self.dds_AOM_A1_Node1.sw.on()
    self.dds_AOM_A1_Node2.sw.on()
    self.dds_AOM_A2_Node1.sw.on()
    self.dds_AOM_A2_Node2.sw.on()
    delay(5 * us)
    self.dds_AOM_A3_Node1.sw.on()
    self.dds_AOM_A3_Node2.sw.on()
    self.dds_AOM_A4_Node1.sw.on()
    self.dds_AOM_A4_Node2.sw.on()
    delay(5 * us)
    ### A5/A6 are the on-chip beams and each node has its own pair, so the
    ### BRANCH is per node too, not just the device name: one node can run with
    ### on-chip beams while the other does not. See point 3 of "three places
    ### two nodes cannot share one line" in the module docstring.
    ### The else branches are REQUIRED, not tidiness. PGC_and_RO_with_on_chip
    ### _beams is True on both nodes, so the `if` never fires and without an
    ### else A5/A6 keep whatever state they were left in -- and the loader
    ### leaves them ON from MOT loading. second_shot DOES have the else and
    ### turns them off, so the first readout ran on six beams and the second on
    ### four. Retention is RO2 vs RO1, so that does not merely add noise, it
    ### makes the ratio meaningless. (The same asymmetry is present in the
    ### single-node original, experiment_functions.py first_shot vs
    ### second_shot; it is NOT introduced by the two-node port.)
    if not self.PGC_and_RO_with_on_chip_beams_Node1:
        self.dds_AOM_A5_Node1.sw.on()
        self.dds_AOM_A6_Node1.sw.on()
    else:
        self.dds_AOM_A5_Node1.sw.off()
        self.dds_AOM_A6_Node1.sw.off()
    if not self.PGC_and_RO_with_on_chip_beams_Node2:
        self.dds_AOM_A5_Node2.sw.on()
        self.dds_AOM_A6_Node2.sw.on()
    else:
        self.dds_AOM_A5_Node2.sw.off()
        self.dds_AOM_A6_Node2.sw.off()
    delay(0.1 * ms)

    with parallel:
        self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_first_shot)
        self.ttl_SPCM1_counter.gate_rising(self.t_SPCM_first_shot)
        self.ttl_SPCM0_OtherNode_counter.gate_rising(self.t_SPCM_first_shot)
        self.ttl_SPCM1_OtherNode_counter.gate_rising(self.t_SPCM_first_shot)

    self.SPCM0_RO1 = self.ttl_SPCM0_counter.fetch_count()
    self.SPCM1_RO1 = self.ttl_SPCM1_counter.fetch_count()
    self.SPCM0_OtherNode_RO1 = self.ttl_SPCM0_OtherNode_counter.fetch_count()
    self.SPCM1_OtherNode_RO1 = self.ttl_SPCM1_OtherNode_counter.fetch_count()
    self.AllSPCMs_RO1 = (self.SPCM0_RO1 + self.SPCM1_RO1
                         + self.SPCM0_OtherNode_RO1 + self.SPCM1_OtherNode_RO1)
    delay(0.1 * ms)
    self.dds_cooling_DP_Node1.sw.off()  ### turn off cooling
    self.dds_cooling_DP_Node2.sw.off()
    self.ttl_repump_switch_Node1.on()  ### turn off MOT RP
    self.ttl_repump_switch_Node2.on()
    delay(5 * us)
    delay(10 * us)


@kernel
def second_shot(self):
    """
    non-chopped second atom readout, both nodes at once.

    Copied from experiment_functions.second_shot. Coils and the two .set()
    calls are per node because every one of those values differs; everything
    else is the original line.

    warning: assumes the fiber AOMs are already on, which is usually the case
    """
    ### set the coils to PGC settings
    self.zotino0_Node1.set_dac(
        [self.AZ_bottom_volts_PGC_Node1, -self.AZ_bottom_volts_PGC_Node1,
         self.AX_volts_PGC_Node1, self.AY_volts_PGC_Node1],
        channels=self.coil_channels_Node1)
    self.zotino0_Node2.set_dac(
        [self.AZ_bottom_volts_PGC_Node2, -self.AZ_bottom_volts_PGC_Node2,
         self.AX_volts_PGC_Node2, self.AY_volts_PGC_Node2],
        channels=self.coil_channels_Node2)
    delay(1 * ms)  ## coils relaxation time

    ### set the FORT AOM to the readout settings
    self.dds_FORT_Node1.set(frequency=self.f_FORT_Node1,
                            amplitude=self.stabilizer_FORT_Node1.amplitudes[1])
    self.dds_FORT_Node2.set(frequency=self.f_FORT_Node2,
                            amplitude=self.stabilizer_FORT_Node2.amplitudes[1])

    ### set the cooling DP AOM to the readout settings
    self.dds_cooling_DP_Node1.set(frequency=self.f_cooling_DP_RO_Node1,
                                  amplitude=self.ampl_cooling_DP_RO_Node1)
    self.dds_cooling_DP_Node2.set(frequency=self.f_cooling_DP_RO_Node2,
                                  amplitude=self.ampl_cooling_DP_RO_Node2)
    delay(5 * us)

    self.ttl_repump_switch_Node1.off()  ### turn on MOT RP
    self.ttl_repump_switch_Node2.off()
    self.dds_cooling_DP_Node1.sw.on()  ### Turn on cooling
    self.dds_cooling_DP_Node2.sw.on()
    delay(5 * us)
    delay(0.1 * ms)

    self.dds_AOM_A1_Node1.sw.on()
    self.dds_AOM_A1_Node2.sw.on()
    self.dds_AOM_A2_Node1.sw.on()
    self.dds_AOM_A2_Node2.sw.on()
    delay(5 * us)
    self.dds_AOM_A3_Node1.sw.on()
    self.dds_AOM_A3_Node2.sw.on()
    self.dds_AOM_A4_Node1.sw.on()
    self.dds_AOM_A4_Node2.sw.on()
    delay(5 * us)
    ### per node, for the same reason as in first_shot
    if not self.PGC_and_RO_with_on_chip_beams_Node1:
        self.dds_AOM_A5_Node1.sw.on()
        self.dds_AOM_A6_Node1.sw.on()
    else:
        self.dds_AOM_A5_Node1.sw.off()
        self.dds_AOM_A6_Node1.sw.off()
        delay(5 * us)
    if not self.PGC_and_RO_with_on_chip_beams_Node2:
        self.dds_AOM_A5_Node2.sw.on()
        self.dds_AOM_A6_Node2.sw.on()
    else:
        self.dds_AOM_A5_Node2.sw.off()
        self.dds_AOM_A6_Node2.sw.off()
        delay(5 * us)
    delay(0.1 * ms)

    with parallel:
        self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_second_shot)
        self.ttl_SPCM1_counter.gate_rising(self.t_SPCM_second_shot)
        self.ttl_SPCM0_OtherNode_counter.gate_rising(self.t_SPCM_second_shot)
        self.ttl_SPCM1_OtherNode_counter.gate_rising(self.t_SPCM_second_shot)

    self.SPCM0_RO2 = self.ttl_SPCM0_counter.fetch_count()
    self.SPCM1_RO2 = self.ttl_SPCM1_counter.fetch_count()
    self.SPCM0_OtherNode_RO2 = self.ttl_SPCM0_OtherNode_counter.fetch_count()
    self.SPCM1_OtherNode_RO2 = self.ttl_SPCM1_OtherNode_counter.fetch_count()
    self.AllSPCMs_RO2 = (self.SPCM0_RO2 + self.SPCM1_RO2
                         + self.SPCM0_OtherNode_RO2 + self.SPCM1_OtherNode_RO2)
    delay(0.1 * ms)
    self.dds_cooling_DP_Node1.sw.off()  ### turn off cooling
    self.dds_cooling_DP_Node2.sw.off()
    self.ttl_repump_switch_Node1.on()  ### turn off MOT RP
    self.ttl_repump_switch_Node2.on()
    delay(5 * us)
    delay(10 * us)


@kernel
def recooling_after_first_shot(self):
    """
    Recool the atoms after the first shot, both nodes.

    THE ONE STAGE THAT CANNOT BE A SINGLE LINE. t_recooling_after_first_shot
    is 0.0 on Node1 and 1.0 ms on Node2, so Node1 does not recool at all while
    Node2 does. Two independent kernels each just delayed its own amount; one
    kernel has to place each node's window from a common t0 and then resync
    past the longer of the two.
    """
    self.dds_cooling_DP_Node1.set(frequency=self.f_cooling_DP_PGC_Node1,
                                  amplitude=self.ampl_cooling_DP_PGC_Node1)
    self.dds_cooling_DP_Node2.set(frequency=self.f_cooling_DP_PGC_Node2,
                                  amplitude=self.ampl_cooling_DP_PGC_Node2)
    delay(0.1 * ms)

    ### Each branch of a `with parallel` starts at the block's entry time and
    ### the cursor afterwards is the MAXIMUM of the branches, which is exactly
    ### "both nodes from a common t0, then resync past the longer". So no
    ### now_mu() anchor and no t_longest bookkeeping -- and no way for that
    ### arithmetic to drift out of step with what the block actually contains.
    with parallel:
        if self.t_recooling_after_first_shot_Node1 > 0.0:
            self.ttl_repump_switch_Node1.off()  ### turn on MOT RP
            self.dds_cooling_DP_Node1.sw.on()
            delay(self.t_recooling_after_first_shot_Node1)
            self.dds_cooling_DP_Node1.sw.off()
            self.ttl_repump_switch_Node1.on()  ### turn off MOT RP

        if self.t_recooling_after_first_shot_Node2 > 0.0:
            self.ttl_repump_switch_Node2.off()
            self.dds_cooling_DP_Node2.sw.on()
            delay(self.t_recooling_after_first_shot_Node2)
            self.dds_cooling_DP_Node2.sw.off()
            self.ttl_repump_switch_Node2.on()

    delay(10 * us)


@kernel
def load_until_atom_in_both_nodes_recycle(self):
    """
    Two-node combined loading.
    * based on load_until_atom_in_both_nodes_recycle, with the TTL
      synchronization removed -- there is one timeline now.
    * Checks for atoms in both nodes at the same time using two_atom_threshold

    Before attempting to load, checks whether atoms are already in the FORTs
    based on RO2. If not, turns on the MOT and FORT light on both nodes and
    watches all four SPCMs. Turns the MOTs off as soon as the summed rate says
    both traps are occupied.
    """

    ### First check if there are already atoms based on RO2
    delay(100 * us)
    atom_loaded = False
    if self.measurement > 0:
        if self.AllSPCMs_RO2/self.t_SPCM_second_shot > self.two_atom_threshold:
            atom_loaded = True

    if not atom_loaded:
        ### Set the coils to MOT loading setting
        self.zotino0_Node1.set_dac(
            [self.AZ_bottom_volts_MOT_Node1, self.AZ_top_volts_MOT_Node1,
             self.AX_volts_MOT_Node1, self.AY_volts_MOT_Node1],
            channels=self.coil_channels_Node1)
        self.zotino0_Node2.set_dac(
            [self.AZ_bottom_volts_MOT_Node2, self.AZ_top_volts_MOT_Node2,
             self.AX_volts_MOT_Node2, self.AY_volts_MOT_Node2],
            channels=self.coil_channels_Node2)
        delay(1 * ms)

        ### set the cooling DP AOM to the MOT settings
        self.dds_cooling_DP_Node1.set(frequency=self.f_cooling_DP_MOT_Node1,
                                      amplitude=self.ampl_cooling_DP_MOT_Node1)
        self.dds_cooling_DP_Node2.set(frequency=self.f_cooling_DP_MOT_Node2,
                                      amplitude=self.ampl_cooling_DP_MOT_Node2)

        self.ttl_repump_switch_Node1.off()  ### turn on MOT RP
        self.ttl_repump_switch_Node2.off()
        delay(5 * us)
        self.dds_cooling_DP_Node1.sw.on()  ### turn on cooling
        self.dds_cooling_DP_Node2.sw.on()

        self.dds_AOM_A1_Node1.sw.on()
        self.dds_AOM_A1_Node2.sw.on()
        delay(5 * us)
        self.dds_AOM_A2_Node1.sw.on()
        self.dds_AOM_A2_Node2.sw.on()
        self.dds_AOM_A3_Node1.sw.on()
        self.dds_AOM_A3_Node2.sw.on()
        delay(5 * us)
        self.dds_AOM_A4_Node1.sw.on()
        self.dds_AOM_A4_Node2.sw.on()
        self.dds_AOM_A5_Node1.sw.on()
        self.dds_AOM_A5_Node2.sw.on()
        delay(5 * us)
        self.dds_AOM_A6_Node1.sw.on()
        self.dds_AOM_A6_Node2.sw.on()
        delay(0.1 * ms)

        ### turn on the FORTs at their loading setpoints
        self.dds_FORT_Node1.set(frequency=self.f_FORT_Node1,
                                amplitude=self.stabilizer_FORT_Node1.amplitude)
        self.dds_FORT_Node2.set(frequency=self.f_FORT_Node2,
                                amplitude=self.stabilizer_FORT_Node2.amplitude)
        self.dds_FORT_Node1.sw.on()
        self.dds_FORT_Node2.sw.on()
        delay(5 * us)

        ### UV trigger is Node1 only; there is no Node2 UV pulse
        self.zotino0_Node1.set_dac([3.5], self.UV_trig_channel_Node1)
        delay(1 * ms)

        t_before_atom = now_mu()
        t_after_atom = now_mu()
        max_tries = self.max_atom_check_tries_two_node
        max_rounds = self.max_loading_rounds_two_node
        atom_check_time = self.t_atom_check_time
        AllSPCMs_atom_check = 0
        AllSPCMs_atom_check_not_loaded = 0
        rounds = 0

        while not atom_loaded and rounds < max_rounds:
            rounds += 1
            ### try_n is reset EVERY round. In the legacy version it was set
            ### once outside the loop and reset only inside the laser-feedback
            ### branch, so with feedback off the inner loop stopped checking
            ### after max_tries and the outer loop spun doing nothing.
            try_n = 0

            while not atom_loaded and try_n < max_tries:
                delay(100 * us)
                with parallel:
                    self.ttl_SPCM0_counter.gate_rising(atom_check_time)
                    self.ttl_SPCM1_counter.gate_rising(atom_check_time)
                    self.ttl_SPCM0_OtherNode_counter.gate_rising(atom_check_time)
                    self.ttl_SPCM1_OtherNode_counter.gate_rising(atom_check_time)

                AllSPCMs_atom_check = (
                    self.ttl_SPCM0_counter.fetch_count()
                    + self.ttl_SPCM1_counter.fetch_count()
                    + self.ttl_SPCM0_OtherNode_counter.fetch_count()
                    + self.ttl_SPCM1_OtherNode_counter.fetch_count())

                try_n += 1
                if try_n == 1:
                    AllSPCMs_atom_check_not_loaded = AllSPCMs_atom_check

                if AllSPCMs_atom_check/atom_check_time > self.two_atom_threshold_for_loading:
                    atom_loaded = True
                    t_after_atom = now_mu()

            if atom_loaded:
                self.set_dataset("time_without_atom", 0.0, broadcast=True)
                self.append_to_dataset("AllSPCMs_atom_check_in_loading", AllSPCMs_atom_check)
                self.append_to_dataset("AllSPCMs_atom_check_in_loading", AllSPCMs_atom_check_not_loaded)
            else:
                #### how long has passed since the previous atom loading
                self.set_dataset(
                    "time_without_atom",
                    self.core.mu_to_seconds(now_mu() - t_before_atom),
                    broadcast=True)

            ### the dataset writes above are host RPCs: they advance wall
            ### clock while the timeline cursor stands still, so without this
            ### the next RTIO event is eventually submitted in the past.
            self.core.break_realtime()

            if not atom_loaded and self.enable_laser_feedback:
                ### reset cooling DP to MOT settings, run feedback, carry on
                self.dds_cooling_DP_Node1.set(frequency=self.f_cooling_DP_MOT_Node1,
                                              amplitude=self.ampl_cooling_DP_MOT_Node1)
                self.dds_cooling_DP_Node2.set(frequency=self.f_cooling_DP_MOT_Node2,
                                              amplitude=self.ampl_cooling_DP_MOT_Node2)
                run_feedback_and_record_FORT_MM_power(self, record_power=False)
                self.n_feedback_per_iteration += 1
                self.dds_FORT_Node1.sw.on()
                self.dds_FORT_Node2.sw.on()
                delay(5 * us)

        if not atom_loaded:
            self.print_async(
                "two-node loading gave up after rounds =", rounds,
                "-- skipping this measurement. Check both MOTs and "
                "two_atom_threshold_for_loading.")

        ### UV off, MOTs away
        self.zotino0_Node1.set_dac([0.0], self.UV_trig_channel_Node1)
        delay(1 * ms)
        self.ttl_repump_switch_Node1.on()  ### turn off MOT RP
        self.ttl_repump_switch_Node2.on()
        self.dds_cooling_DP_Node1.sw.off()  ### turn off cooling
        self.dds_cooling_DP_Node2.sw.off()
        delay(5 * us)
        delay(self.t_MOT_dissipation)

        ### set the coils to PGC settings; effectively turns the coils off
        self.zotino0_Node1.set_dac(
            [self.AZ_bottom_volts_PGC_Node1, -self.AZ_bottom_volts_PGC_Node1,
             self.AX_volts_PGC_Node1, self.AY_volts_PGC_Node1],
            channels=self.coil_channels_Node1)
        self.zotino0_Node2.set_dac(
            [self.AZ_bottom_volts_PGC_Node2, -self.AZ_bottom_volts_PGC_Node2,
             self.AX_volts_PGC_Node2, self.AY_volts_PGC_Node2],
            channels=self.coil_channels_Node2)
        delay(1 * ms)

        self.atom_loading_time = self.core.mu_to_seconds(t_after_atom - t_before_atom)
        self.append_to_dataset("Atom_loading_time", self.atom_loading_time)
        self.append_to_dataset("atom_loading_wall_clock", now_mu())
        if atom_loaded:
            self.n_atom_loaded_per_iteration += 1
        delay(1 * ms)
        self.core.break_realtime()

    ###########  PGC on the trapped atoms  #############
    ### do_PGC_after_loading is per node, not shared: each node's PGC window is
    ### its own stage, and t_PGC_after_loading already differs between them
    ### (0.6 vs 1.0 ms). A node that is not doing PGC simply has no window.
    if self.do_PGC_after_loading_Node1 or self.do_PGC_after_loading_Node2:
        ### THE OTHER STAGE THAT CANNOT BE A SINGLE LINE: t_PGC_after_loading
        ### is 0.6 ms on Node1 and 1.0 ms on Node2, so each node's PGC window
        ### is placed from a common t0 and the cursor resynced past the longer.
        self.dds_cooling_DP_Node1.set(frequency=self.f_cooling_DP_PGC_Node1,
                                      amplitude=self.ampl_cooling_DP_PGC_Node1)
        self.dds_cooling_DP_Node2.set(frequency=self.f_cooling_DP_PGC_Node2,
                                      amplitude=self.ampl_cooling_DP_PGC_Node2)

        self.dds_AOM_A1_Node1.sw.on()
        self.dds_AOM_A1_Node2.sw.on()
        delay(5 * us)
        self.dds_AOM_A2_Node1.sw.on()
        self.dds_AOM_A2_Node2.sw.on()
        self.dds_AOM_A3_Node1.sw.on()
        self.dds_AOM_A3_Node2.sw.on()
        delay(5 * us)
        self.dds_AOM_A4_Node1.sw.on()
        self.dds_AOM_A4_Node2.sw.on()
        delay(0.1 * ms)
        ### on-chip beams: each node's own pair, on its own flag
        if (self.do_PGC_after_loading_Node1
                and not self.PGC_and_RO_with_on_chip_beams_Node1):
            self.dds_AOM_A5_Node1.sw.on()
            self.dds_AOM_A6_Node1.sw.on()
        if (self.do_PGC_after_loading_Node2
                and not self.PGC_and_RO_with_on_chip_beams_Node2):
            self.dds_AOM_A5_Node2.sw.on()
            self.dds_AOM_A6_Node2.sw.on()
            delay(5 * us)

        ### `with parallel` does BOTH jobs the at_mu bookkeeping used to do:
        ### every statement in the block starts at the block's entry time, and
        ### the cursor afterwards is the MAXIMUM of the branches' end times. So
        ### the common t0 and the resync-past-the-longer are both automatic and
        ### t_pgc_longest is gone -- along with the bug class where that
        ### arithmetic drifts out of step with the block's contents.
        ###
        ### Note it is a timeline construct, not concurrency: the CPU still
        ### submits both branches' events one after another, and both branches'
        ### FIRST events land on the same timestamp. Two events here, so well
        ### inside the four-per-timestamp budget.
        with parallel:
            if self.do_PGC_after_loading_Node1:
                self.ttl_repump_switch_Node1.off()  ### turn on MOT RP
                self.dds_cooling_DP_Node1.sw.on()
                delay(self.t_PGC_after_loading_Node1)
                self.dds_cooling_DP_Node1.sw.off()
                self.ttl_repump_switch_Node1.on()

            if self.do_PGC_after_loading_Node2:
                self.ttl_repump_switch_Node2.off()
                self.dds_cooling_DP_Node2.sw.on()
                delay(self.t_PGC_after_loading_Node2)
                self.dds_cooling_DP_Node2.sw.off()
                self.ttl_repump_switch_Node2.on()

        delay(10 * us)


@kernel
def two_node_alternating_shot(self):
    """
    Alternating two-node readout: even windows Node1, odd windows Node2.

    Still needed, and unchanged in purpose by the gateware: the beamsplitter
    fan-out puts both nodes' fluorescence on all four SPCMs, so separating
    them in time is the only way to attribute a photon to a trap. All four
    counters gate in every window.

    One kernel owns the window parity now. In the legacy version both cores
    ran this whole loop and each only gated its own light, so both nodes
    accumulated into both _alice and _bob and the attribution depended on two
    independent cores agreeing on parity.

    Windows are placed absolutely from one t0, and every fetch_count happens
    after the loop. The legacy version blocked on four fetch_count calls per
    window, which drops the loop to ~zero slack from window 2 onward with
    Node2 DDS writes still inside it.

    LIMIT ON n_alternating_RO_windows_per_node: deferring the fetches means
    every window's count sits in the RTIO input FIFO until the loop ends, and
    each counter queues one total per window -- so 2 * n totals per counter,
    all four counters in parallel. At the default n = 10 that is 20 per
    counter, comfortably inside the 64-deep Kasli input FIFO. Past n = 32 it
    overflows. The failure is loud (RTIOOverflow, not a wrong number), and the
    fix is to fetch in chunks rather than to re-block per window.
    """
    self.zotino0_Node1.set_dac(
        [self.AZ_bottom_volts_PGC_Node1, -self.AZ_bottom_volts_PGC_Node1,
         self.AX_volts_PGC_Node1, self.AY_volts_PGC_Node1],
        channels=self.coil_channels_Node1)
    self.zotino0_Node2.set_dac(
        [self.AZ_bottom_volts_PGC_Node2, -self.AZ_bottom_volts_PGC_Node2,
         self.AX_volts_PGC_Node2, self.AY_volts_PGC_Node2],
        channels=self.coil_channels_Node2)
    delay(1 * ms)

    self.dds_FORT_Node1.set(frequency=self.f_FORT_Node1,
                            amplitude=self.stabilizer_FORT_Node1.amplitudes[1])
    self.dds_FORT_Node2.set(frequency=self.f_FORT_Node2,
                            amplitude=self.stabilizer_FORT_Node2.amplitudes[1])
    self.dds_cooling_DP_Node1.set(frequency=self.f_cooling_DP_RO_Node1,
                                  amplitude=self.ampl_cooling_DP_RO_Node1)
    self.dds_cooling_DP_Node2.set(frequency=self.f_cooling_DP_RO_Node2,
                                  amplitude=self.ampl_cooling_DP_RO_Node2)
    delay(5 * us)

    self.dds_AOM_A1_Node1.sw.on()
    self.dds_AOM_A1_Node2.sw.on()
    self.dds_AOM_A2_Node1.sw.on()
    self.dds_AOM_A2_Node2.sw.on()
    delay(5 * us)
    self.dds_AOM_A3_Node1.sw.on()
    self.dds_AOM_A3_Node2.sw.on()
    self.dds_AOM_A4_Node1.sw.on()
    self.dds_AOM_A4_Node2.sw.on()
    delay(5 * us)
    ### per node, for the same reason as in first_shot
    if not self.PGC_and_RO_with_on_chip_beams_Node1:
        self.dds_AOM_A5_Node1.sw.on()
        self.dds_AOM_A6_Node1.sw.on()
    else:
        self.dds_AOM_A5_Node1.sw.off()
        self.dds_AOM_A6_Node1.sw.off()
        delay(5 * us)
    if not self.PGC_and_RO_with_on_chip_beams_Node2:
        self.dds_AOM_A5_Node2.sw.on()
        self.dds_AOM_A6_Node2.sw.on()
    else:
        self.dds_AOM_A5_Node2.sw.off()
        self.dds_AOM_A6_Node2.sw.off()
        delay(5 * us)
    delay(0.1 * ms)

    ### start dark
    self.ttl_repump_switch_Node1.on()
    self.ttl_repump_switch_Node2.on()
    self.dds_cooling_DP_Node1.sw.off()
    self.dds_cooling_DP_Node2.sw.off()
    delay(5 * us)

    n_windows = self.n_alternating_RO_windows_per_node
    pad_mu = self.core.seconds_to_mu(self.t_alternating_RO_pad)
    window_mu = self.core.seconds_to_mu(self.t_alternating_RO_window)
    node2_offset_mu = self.t_Node2_rtio_offset_mu
    period_mu = window_mu + 3 * pad_mu

    self.core.break_realtime()
    t0 = now_mu()

    for i in range(2 * n_windows):
        t_window = t0 + i * period_mu
        alice_window = (i % 2) == 0

        ### only the active node turns on readout light
        if alice_window:
            at_mu(t_window)
            self.ttl_repump_switch_Node1.off()
            self.dds_cooling_DP_Node1.sw.on()
        else:
            at_mu(t_window - node2_offset_mu)
            self.ttl_repump_switch_Node2.off()
            self.dds_cooling_DP_Node2.sw.on()

        ### the counters are master-local, so never offset
        at_mu(t_window + pad_mu)
        with parallel:
            self.ttl_SPCM0_counter.gate_rising_mu(window_mu)
            self.ttl_SPCM1_counter.gate_rising_mu(window_mu)
            self.ttl_SPCM0_OtherNode_counter.gate_rising_mu(window_mu)
            self.ttl_SPCM1_OtherNode_counter.gate_rising_mu(window_mu)

        if alice_window:
            at_mu(t_window + pad_mu + window_mu)
            self.dds_cooling_DP_Node1.sw.off()
            self.ttl_repump_switch_Node1.on()
        else:
            at_mu(t_window + pad_mu + window_mu - node2_offset_mu)
            self.dds_cooling_DP_Node2.sw.off()
            self.ttl_repump_switch_Node2.on()

    at_mu(t0 + 2 * n_windows * period_mu)
    delay(1 * ms)

    ### deferred fetches, in emission order
    self.AllSPCMs_alternating_RO_alice = 0
    self.AllSPCMs_alternating_RO_bob = 0
    for i in range(2 * n_windows):
        window_count = (self.ttl_SPCM0_counter.fetch_count()
                        + self.ttl_SPCM1_counter.fetch_count()
                        + self.ttl_SPCM0_OtherNode_counter.fetch_count()
                        + self.ttl_SPCM1_OtherNode_counter.fetch_count())
        if (i % 2) == 0:
            self.AllSPCMs_alternating_RO_alice += window_count
        else:
            self.AllSPCMs_alternating_RO_bob += window_count
    self.core.break_realtime()


@kernel
def end_measurement(self):
    """
    End the measurement by setting datasets and deciding whether to increment
    the measurement index.

    Copied from experiment_functions.end_measurement, with ONE difference: the
    magnetometer block is gone. It reads

        if self.monitor_magnetometer_in_end_measurement:
            if self.which_node == "alice": measure_Magnetometer(self)
            else:                          measure_Magnetometer_Node2_X_Y_Z(self)

    and although that is dead in a two-node run (the flag is False on both
    nodes, and a two-node sequence samples no magnetometers), the ARTIQ
    compiler types every branch it can reach -- so it still demands bare
    coil_channels, zotino0 and Magnetometer_*_ch, and coil_channels genuinely
    differs per node ([0,1,13,14] vs [0,1,2,3]). That one dead branch is why
    this is a copy rather than an import.

    Everything else is the original line, including the dataset names, so the
    applets and the analysis notebooks read a two-node run exactly as they
    read a single-node one.
    """
    in_health_check = self.in_health_check

    ### update the datasets
    self.set_dataset(self.measurements_progress, 100 * self.measurement / self.n_measurements, broadcast=True)

    self.append_to_dataset('SPCM0_RO1_current_iteration', self.SPCM0_RO1)
    self.append_to_dataset('SPCM1_RO1_current_iteration', self.SPCM1_RO1)
    self.append_to_dataset('SPCM0_OtherNode_RO1_current_iteration', self.SPCM0_OtherNode_RO1)
    self.append_to_dataset('SPCM1_OtherNode_RO1_current_iteration', self.SPCM1_OtherNode_RO1)
    delay(1 * ms)
    self.append_to_dataset('SPCM0_RO2_current_iteration', self.SPCM0_RO2)
    self.append_to_dataset('SPCM1_RO2_current_iteration', self.SPCM1_RO2)
    self.append_to_dataset('SPCM0_OtherNode_RO2_current_iteration', self.SPCM0_OtherNode_RO2)
    self.append_to_dataset('SPCM1_OtherNode_RO2_current_iteration', self.SPCM1_OtherNode_RO2)
    delay(1 * ms)
    self.append_to_dataset('AllSPCMs_RO1_current_iteration', self.AllSPCMs_RO1)
    self.append_to_dataset('AllSPCMs_RO2_current_iteration', self.AllSPCMs_RO2)
    delay(1 * ms)
    ### Alternating RO
    self.append_to_dataset('AllSPCMs_alternating_RO_alice_current_iteration', self.AllSPCMs_alternating_RO_alice)
    self.append_to_dataset('AllSPCMs_alternating_RO_bob_current_iteration', self.AllSPCMs_alternating_RO_bob)
    delay(1*ms)

    self.SPCM0_RO1_list[self.measurement] = self.SPCM0_RO1
    self.SPCM1_RO1_list[self.measurement] = self.SPCM1_RO1
    self.SPCM0_RO2_list[self.measurement] = self.SPCM0_RO2
    self.SPCM1_RO2_list[self.measurement] = self.SPCM1_RO2

    self.SPCM0_OtherNode_RO1_list[self.measurement] = self.SPCM0_OtherNode_RO1
    self.SPCM1_OtherNode_RO1_list[self.measurement] = self.SPCM1_OtherNode_RO1
    self.SPCM0_OtherNode_RO2_list[self.measurement] = self.SPCM0_OtherNode_RO2
    self.SPCM1_OtherNode_RO2_list[self.measurement] = self.SPCM1_OtherNode_RO2
    delay(1 * ms)
    self.AllSPCMs_RO1_list[self.measurement] = self.AllSPCMs_RO1
    self.AllSPCMs_RO2_list[self.measurement] = self.AllSPCMs_RO2
    self.atom_loading_time_list[self.measurement] = self.atom_loading_time

    delay(1*ms)

    self.append_to_dataset('AllSPCMs_alternating_RO_alice', self.AllSPCMs_alternating_RO_alice)
    self.append_to_dataset('AllSPCMs_alternating_RO_bob', self.AllSPCMs_alternating_RO_bob)
    delay(1*ms)

    self.measurement += 1
    delay(1 * ms)
    if not in_health_check:  ## advance and in_health_check are different type so can't be mixed.
        self.append_to_dataset('SPCM0_RO1', self.SPCM0_RO1)
        self.append_to_dataset('SPCM1_RO1', self.SPCM1_RO1)
        self.append_to_dataset('SPCM0_RO2', self.SPCM0_RO2)
        self.append_to_dataset('SPCM1_RO2', self.SPCM1_RO2)
        delay(1 * ms)
        self.append_to_dataset('SPCM0_OtherNode_RO1', self.SPCM0_OtherNode_RO1)
        self.append_to_dataset('SPCM1_OtherNode_RO1', self.SPCM1_OtherNode_RO1)
        self.append_to_dataset('SPCM0_OtherNode_RO2', self.SPCM0_OtherNode_RO2)
        self.append_to_dataset('SPCM1_OtherNode_RO2', self.SPCM1_OtherNode_RO2)

        self.append_to_dataset('AllSPCMs_RO1', self.AllSPCMs_RO1)
        self.append_to_dataset('AllSPCMs_RO2', self.AllSPCMs_RO2)
        delay(1 * ms)
    else:
        self.append_to_dataset('SPCM0_RO1_in_health_check', self.SPCM0_RO1)
        self.append_to_dataset('SPCM1_RO1_in_health_check', self.SPCM1_RO1)
        self.append_to_dataset('SPCM0_RO2_in_health_check', self.SPCM0_RO2)
        self.append_to_dataset('SPCM1_RO2_in_health_check', self.SPCM1_RO2)
        delay(1 * ms)
        self.append_to_dataset('SPCM0_OtherNode_RO1_in_health_check', self.SPCM0_OtherNode_RO1)
        self.append_to_dataset('SPCM1_OtherNode_RO1_in_health_check', self.SPCM1_OtherNode_RO1)
        self.append_to_dataset('SPCM0_OtherNode_RO2_in_health_check', self.SPCM0_OtherNode_RO2)
        self.append_to_dataset('SPCM1_OtherNode_RO2_in_health_check', self.SPCM1_OtherNode_RO2)

        self.append_to_dataset('AllSPCMs_RO1_in_health_check', self.AllSPCMs_RO1)
        self.append_to_dataset('AllSPCMs_RO2_in_health_check', self.AllSPCMs_RO2)

    ### ~25 host RPCs against ~10 ms of coded delay leaves the cursor behind
    ### the wall clock, and the next RTIO event -- in the CALLER -- underflows.
    self.core.break_realtime()


@kernel
def Two_nodes_atom_loading_experiment(self):
    """
    ** modified for master-satellite **
    Simple atom loading experiment in both nodes
    - checking for atoms in both nodes at the same time
    - based on atom_loading_2_experiment

    Sequence as of 2026.07.13
    1. Load atoms in both nodes simultaneously
    2. First shot
    3. FORT drop if > 0
    4. testing individual node shot - Alternating RO
    5. Second shot
    6. end_measurement
    """

    self.core.reset()
    self.require_D1_lock_to_advance = False  # override experiment variable

    self.n_feedback_per_iteration = 2  ### number of times the feedback runs in each iteration. Updates in atom loading subroutines.
    ### Required only for averaging RF powers over iterations in analysis. Starts with 2 because RF is measured at least 2 times
    ### in each iteration.
    self.n_atom_loaded_per_iteration = 0

    if self.enable_laser_feedback:
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP_Node1.set(frequency=self.f_cooling_DP_MOT_Node1,
                                      amplitude=self.ampl_cooling_DP_MOT_Node1)
        self.dds_cooling_DP_Node2.set(frequency=self.f_cooling_DP_MOT_Node2,
                                      amplitude=self.ampl_cooling_DP_MOT_Node2)
        run_feedback_and_record_FORT_MM_power(self)

    self.measurement = 0
    while self.measurement < self.n_measurements:
        delay(10 * ms)

        load_until_atom_in_both_nodes_recycle(self)

        delay(1 * ms)

        first_shot(self)
        delay(1 * ms)

        if self.t_recooling_after_first_shot_Node1 > 0 or self.t_recooling_after_first_shot_Node2 > 0:
            recooling_after_first_shot(self)

        ### first_shot doesn't turn off the fiber AOMs. thus, PR was actually being done with all 6 beams!!!! :(
        self.dds_AOM_A1_Node1.sw.off()
        self.dds_AOM_A1_Node2.sw.off()
        self.dds_AOM_A2_Node1.sw.off()
        self.dds_AOM_A2_Node2.sw.off()
        delay(5 * us)
        self.dds_AOM_A3_Node1.sw.off()
        self.dds_AOM_A3_Node2.sw.off()
        self.dds_AOM_A4_Node1.sw.off()
        self.dds_AOM_A4_Node2.sw.off()
        delay(5 * us)
        delay(0.1 * ms)
        ### per node, for the same reason as in first_shot
        if not self.PGC_and_RO_with_on_chip_beams_Node1:
            self.dds_AOM_A5_Node1.sw.off()
            self.dds_AOM_A6_Node1.sw.off()
        if not self.PGC_and_RO_with_on_chip_beams_Node2:
            self.dds_AOM_A5_Node2.sw.off()
            self.dds_AOM_A6_Node2.sw.off()
            delay(5 * us)

        ### The FORT drop is the knob retention is measured against, and
        ### t_FORT_drop is per node: the two drops are placed from a common t0
        ### and the cursor resynced past the longer by `with parallel`, so a
        ### node with a 0 drop simply has no gap and each node still gets
        ### exactly its own t_FORT_drop. Sharing one bare value would hand Node2
        ### Node1's drop and quietly make the two nodes' retention numbers
        ### incomparable -- the one thing a two-node run exists to compare.
        with parallel:
            if self.t_FORT_drop_Node1 > 0.0:
                self.dds_FORT_Node1.sw.off()
                delay(self.t_FORT_drop_Node1)
                self.dds_FORT_Node1.sw.on()

            if self.t_FORT_drop_Node2 > 0.0:
                self.dds_FORT_Node2.sw.off()
                delay(self.t_FORT_drop_Node2)
                self.dds_FORT_Node2.sw.on()

        delay(self.t_delay_between_shots)

        # delay(1*ms)
        # two_node_alternating_shot(self)
        # delay(1 * ms)

        second_shot(self)

        end_measurement(self)

    self.append_to_dataset('n_feedback_per_iteration', self.n_feedback_per_iteration)
    self.append_to_dataset('n_atom_loaded_per_iteration', self.n_atom_loaded_per_iteration)


@kernel
def Two_nodes_alternating_shot_experiment(self):
    """
    ** master-satellite **
    Load atoms in both nodes, then attribute photons to each node by taking
    the readout windows alternately instead of jointly.

    This is Two_nodes_atom_loading_experiment with two_node_alternating_shot in
    place of the joint second shot. It exists as its own registered experiment
    rather than as a flag inside the loading sequence for two reasons:

      * a flag would have to arrive through the override dictionary, and a
        missing flag fails to COMPILE, not to run;
      * two_node_alternating_shot is otherwise called by nothing, and ARTIQ
        only type-checks functions reachable from the kernel entry point -- so
        without this harness those ~110 lines of kernel code are never
        compiled by anything, and a compile sweep that "passes" says nothing
        about them.

    The joint first shot still runs, so AllSPCMs_RO1 is directly comparable
    with a joint run; it is only the attribution step that changes.
    """

    self.core.reset()
    self.require_D1_lock_to_advance = False  # override experiment variable

    self.n_feedback_per_iteration = 2
    self.n_atom_loaded_per_iteration = 0

    if self.enable_laser_feedback:
        self.dds_cooling_DP_Node1.set(frequency=self.f_cooling_DP_MOT_Node1,
                                      amplitude=self.ampl_cooling_DP_MOT_Node1)
        self.dds_cooling_DP_Node2.set(frequency=self.f_cooling_DP_MOT_Node2,
                                      amplitude=self.ampl_cooling_DP_MOT_Node2)
        run_feedback_and_record_FORT_MM_power(self)

    self.measurement = 0
    while self.measurement < self.n_measurements:
        delay(10 * ms)

        load_until_atom_in_both_nodes_recycle(self)

        delay(1 * ms)

        first_shot(self)
        delay(1 * ms)

        ### first_shot leaves the fiber AOMs on
        self.dds_AOM_A1_Node1.sw.off()
        self.dds_AOM_A1_Node2.sw.off()
        self.dds_AOM_A2_Node1.sw.off()
        self.dds_AOM_A2_Node2.sw.off()
        delay(5 * us)
        self.dds_AOM_A3_Node1.sw.off()
        self.dds_AOM_A3_Node2.sw.off()
        self.dds_AOM_A4_Node1.sw.off()
        self.dds_AOM_A4_Node2.sw.off()
        delay(5 * us)
        delay(0.1 * ms)
        if not self.PGC_and_RO_with_on_chip_beams_Node1:
            self.dds_AOM_A5_Node1.sw.off()
            self.dds_AOM_A6_Node1.sw.off()
        if not self.PGC_and_RO_with_on_chip_beams_Node2:
            self.dds_AOM_A5_Node2.sw.off()
            self.dds_AOM_A6_Node2.sw.off()
            delay(5 * us)

        delay(self.t_delay_between_shots)

        ### the attribution step, in place of a joint second shot
        two_node_alternating_shot(self)
        delay(1 * ms)

        ### second_shot still runs so the retention criterion and every
        ### downstream applet and analysis script see the dataset they expect;
        ### the alternating counts are recorded alongside it.
        second_shot(self)

        ### end_measurement already appends AllSPCMs_alternating_RO_alice and
        ### _bob, and their _current_iteration variants, so do NOT append them
        ### here: that would put two entries per measurement into datasets the
        ### applets and analysis zip against the one-per-measurement ones.
        ### It also means a joint run and an alternating run produce
        ### identically shaped h5 files -- the joint one just leaves these two
        ### at zero.
        end_measurement(self)

    self.append_to_dataset('n_feedback_per_iteration', self.n_feedback_per_iteration)
    self.append_to_dataset('n_atom_loaded_per_iteration', self.n_atom_loaded_per_iteration)
