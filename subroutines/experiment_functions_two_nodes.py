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

SHARED SETTINGS are read bare, and that is about values rather than devices.
Two kinds, both declared in ExperimentVariables_master_satellite_global.py:
  * one value for both nodes in every mode (t_delay_between_shots,
    t_MOT_dissipation, enable_laser_feedback): one interval on the one shared
    timeline, or a run-wide switch, so there is only one value to have;
  * the two-node sequence settings (the joint gates t_SPCM_first_shot,
    t_SPCM_second_shot and t_atom_check_time -- ONE gate of all four
    master-local SPCMs each -- and n_excitation_attempts, t_pumping,
    t_photon_collection_time, ...): Base publishes self.<x> from the global
    <x>_two_nodes in two_nodes mode, leaving the per-node values single-node
    experiments use alone. See BaseExperimentMasterSatellite.TWO_NODES_SUFFIX.

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

### The kernel version, as in experiment_functions.py. Needed since the first
### recording that sets a DDS from a dBm value (record_CW_optical_pumping_both
### _nodes) -- which nothing called until Two_nodes_HOM_experiment, so the
### missing import never surfaced.
from utilities.conversions import dB_to_V_kernel as dB_to_V


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
def record_CW_optical_pumping_both_nodes(self):
    """Record ONE DMA trace that optically pumps both nodes together.

    Ported from record_CW_optical_pumping_node1 (now in
    experiment_functions_obsolete.py) and record_CW_optical_pumping_node2 (now
    record_CW_optical_pumping in experiment_functions.py), merged rather than
    replayed side by side.

    WHY ONE TRACE AND NOT TWO PLAYED IN PARALLEL. core_dma.playback_handle
    ends with delay_mu(advance_mu), so two handles played sequentially run back
    to back, and inside a `with parallel` the DMA engine submits ALL of trace 1
    and then restarts trace 2 at the block origin -- one rewind spanning an
    entire trace, which the SED has to absorb across its 8 lanes
    (gateware/rtio/sed/__init__.py:30 and ":50"). One trace whose timestamps
    only increase has no rewind at all, and every stage boundary becomes a
    single timestamp that both nodes share. The DMA reader is a plain
    sequential memory stream with no sorting (gateware/rtio/dma.py), so record
    order IS playback order; a `with parallel` INSIDE the recording would
    reintroduce the same rewind, which is why there is none below. (The lanes
    are per crate: the master's SED sees only Node1's events and the
    satellite's SED only Node2's, so a rewind costs lanes only on a crate whose
    OWN events go backwards -- gateware/rtio/cri.py routes each write to one
    destination.)

    It reaches both crates because the master's DMA CRI is fanned out to the
    local RTIO core and every DRTIO destination
    (gateware/targets/kasli.py:394-396), and DRTIO costs nothing per event
    (measured 2026-09-10, RID 38429: 713 mu master vs 711 mu satellite).

    WHAT "PARALLEL" MEANS HERE. Both nodes' events for the same step are
    emitted at the SAME timestamp, in groups of at most four, with 5 us
    between groups. Four per timestamp is the house rule, well inside the 8
    SED lanes beyond which events are silently discarded. The pumping window
    opens and closes on one timestamp for both nodes, and both Zotinos latch
    on one timestamp, so the stages really are simultaneous.

    STAGE DURATIONS ARE SHARED: t_pumping, t_depumping and t_OP_coil_relax_time
    are the two-node settings (t_pumping_two_nodes etc., read here under the
    bare names), one value for both nodes, since both nodes' stages open and
    close on the same timestamps.

    NO t_Node2_rtio_offset_mu CORRECTION, deliberately. That cancels the
    electrical skew between the two crates' output pins and only matters where
    a window is microseconds wide -- the alternating readout, which validates
    |offset| < t_alternating_RO_pad for exactly that reason. Every stage here
    is 30 us or longer, and this file's loading and PGC sequences do not apply
    it either.

    THE NODES ARE WIRED ALIKE since the 2026-10-20 Node1 rewire: on each, D1
    comes from an external source gated by ttl_D1_pumping_switch_NodeX and
    reaches the atom through the GRIN1 AOM, so dds_GRIN1_NodeX runs at the D1
    pumping tone and ttl_GRIN1_switch_NodeX opens with the D1 gate. Both nodes
    deliver the pumping repump through fiber AOMs A5/A6 only: Node1's A1-A4
    were dropped here on 2026-10-23, following the single-node
    record_CW_optical_pumping after the user's 2026-10-21 test on Node1.
    """
    ### The shared two-node settings (t_OP_coil_relax_time_two_nodes etc.)
    t_coil_relax = self.t_OP_coil_relax_time
    t_pumping = self.t_pumping
    t_depumping = self.t_depumping

    with self.core_dma.record("CW_optical_pumping_both_nodes"):
        ### Turn off unused light. Eight events, so two groups.
        self.ttl_repump_switch_Node1.on()
        self.ttl_repump_switch_Node2.on()
        self.ttl_exc0_switch_Node1.on()
        self.ttl_exc0_switch_Node2.on()
        delay(5 * us)
        self.dds_cooling_DP_Node1.sw.off()
        self.dds_cooling_DP_Node2.sw.off()
        ### make sure both dedicated D1 gates start closed
        self.ttl_D1_pumping_switch_Node1.on()
        self.ttl_D1_pumping_switch_Node2.on()
        delay(5 * us)

        ### Configure the pumping-repump delivery AOMs. The four-per-group rule
        ### covers set() as well as the switches -- within a node these share a
        ### Urukul SPI bus, and the lab rule names dds.set explicitly.
        self.dds_AOM_A5_Node1.set(frequency=self.AOM_A5_freq_Node1, amplitude=dB_to_V(-5.0))
        self.dds_AOM_A6_Node1.set(frequency=self.AOM_A6_freq_Node1, amplitude=dB_to_V(-5.0))
        self.dds_AOM_A5_Node2.set(frequency=self.AOM_A5_freq_Node2, amplitude=dB_to_V(-5.0))
        self.dds_AOM_A6_Node2.set(frequency=self.AOM_A6_freq_Node2, amplitude=dB_to_V(-5.0))
        delay(5 * us)

        ### ... and switch them on. A5/A6 fire on both nodes at one timestamp.
        self.dds_AOM_A5_Node1.sw.on()
        self.dds_AOM_A6_Node1.sw.on()
        self.dds_AOM_A5_Node2.sw.on()
        self.dds_AOM_A6_Node2.sw.on()
        delay(5 * us)

        ### Program both nodes' coil DACs forward in time, then LOAD both at
        ### ONE timestamp. load() is the atomic update, so the shift-register
        ### writes stay sequential while the coils still move together -- and
        ### these two load() pulses are the thing to put on a scope to prove
        ### this function does what it claims. Separate Zotinos on separate
        ### crates, so the eight writes contend for nothing.
        self.zotino0_Node1.write_dac(self.coil_channels_Node1[0], self.AZ_bottom_volts_OP_Node1)
        self.zotino0_Node1.write_dac(self.coil_channels_Node1[1], -self.AZ_bottom_volts_OP_Node1)
        self.zotino0_Node1.write_dac(self.coil_channels_Node1[2], self.AX_volts_OP_Node1)
        self.zotino0_Node1.write_dac(self.coil_channels_Node1[3], self.AY_volts_OP_Node1)
        self.zotino0_Node2.write_dac(self.coil_channels_Node2[0], self.AZ_bottom_volts_OP_Node2)
        self.zotino0_Node2.write_dac(self.coil_channels_Node2[1], -self.AZ_bottom_volts_OP_Node2)
        self.zotino0_Node2.write_dac(self.coil_channels_Node2[2], self.AX_volts_OP_Node2)
        self.zotino0_Node2.write_dac(self.coil_channels_Node2[3], self.AY_volts_OP_Node2)
        delay(1.5 * us)
        self.zotino0_Node1.load()
        self.zotino0_Node2.load()

        ### Physical coil-relaxation time, shared by both nodes
        delay(t_coil_relax)

        ### Both nodes' D1 reaches the atom through GRIN1, so the GRIN1 tone
        ### must be up before the gates open.
        self.dds_GRIN1_Node1.set(frequency=self.f_GRIN1_D1_pumping_Node1,
                                 amplitude=dB_to_V(self.p_GRIN1_D1_pumping_Node1))
        self.dds_GRIN1_Node2.set(frequency=self.f_GRIN1_D1_pumping_Node2,
                                 amplitude=dB_to_V(self.p_GRIN1_D1_pumping_Node2))
        self.dds_GRIN1_Node1.sw.on()
        self.dds_GRIN1_Node2.sw.on()
        delay(2 * us)

        ### Pumping repump on, 5 us AHEAD of the D1 light rather than with it,
        ### so that the pumping start below gets a timestamp to itself. The
        ### repump only has to be on; it does not define the window.
        self.ttl_pumping_repump_switch_Node1.off()
        self.ttl_pumping_repump_switch_Node2.off()
        delay(5 * us)

        ### OPTICAL PUMPING STARTS HERE -- one timestamp, both nodes
        self.ttl_D1_pumping_switch_Node1.off()
        self.ttl_D1_pumping_switch_Node2.off()
        self.ttl_GRIN1_switch_Node1.off()
        self.ttl_GRIN1_switch_Node2.off()

        delay(t_pumping)

        ### ... AND ENDS HERE -- one timestamp, both nodes. The D1 gates close
        ### and GRIN1's RF goes off together, as in the single-node
        ### record_CW_optical_pumping since e88bc0a (2026-10-26): before that
        ### the GRIN1 RF stayed on behind its closed gate through depumping.
        ### The GRIN1 gates close t_depumping later instead, with the pumping
        ### repump -- by then there is no RF behind them -- which keeps this
        ### timestamp at four events and the trace the same length.
        self.ttl_D1_pumping_switch_Node1.on()
        self.ttl_D1_pumping_switch_Node2.on()
        self.dds_GRIN1_Node1.sw.off()
        self.dds_GRIN1_Node2.sw.off()

        delay(t_depumping)

        ### Turn off every pumping RF and light path. Eight events, so two
        ### groups. The pumping repump stays on through depumping, as in both
        ### single-node originals, and goes off here.
        self.ttl_pumping_repump_switch_Node1.on()
        self.ttl_pumping_repump_switch_Node2.on()
        self.ttl_GRIN1_switch_Node1.on()
        self.ttl_GRIN1_switch_Node2.on()
        delay(5 * us)
        self.dds_AOM_A5_Node1.sw.off()
        self.dds_AOM_A6_Node1.sw.off()
        self.dds_AOM_A5_Node2.sw.off()
        self.dds_AOM_A6_Node2.sw.off()
        delay(5 * us)

        ### Restore the feedback amplitudes while the RF switches are off.
        ### .amplitude is amplitudes[0], kept in step by the feedback routine
        ### (aom_feedback.py:182) -- the loading setpoint, as in the originals.
        self.dds_AOM_A5_Node1.set(frequency=self.AOM_A5_freq_Node1, amplitude=self.stabilizer_AOM_A5_Node1.amplitude)
        self.dds_AOM_A6_Node1.set(frequency=self.AOM_A6_freq_Node1, amplitude=self.stabilizer_AOM_A6_Node1.amplitude)
        self.dds_AOM_A5_Node2.set(frequency=self.AOM_A5_freq_Node2, amplitude=self.stabilizer_AOM_A5_Node2.amplitude)
        self.dds_AOM_A6_Node2.set(frequency=self.AOM_A6_freq_Node2, amplitude=self.stabilizer_AOM_A6_Node2.amplitude)


@kernel
def record_recooling_both_nodes(self):
    """Record ONE DMA trace "recooling_both_nodes" that recools both atoms.

    Two-node version of record_recooling (experiment_functions.py), merged the
    way record_CW_optical_pumping_both_nodes is: one trace whose timestamps
    only increase, both nodes' events for a step on one timestamp, at most four
    events per timestamp with 5 us between groups.

    The single-node original switches eight events on ONE timestamp (cooling,
    repump, A1-A4, A5, A6) -- the SED's hard limit -- so a doubled copy would
    drop events. Here the fiber AOMs are switched first, in groups, and the
    light then reaches both atoms on ONE timestamp, when the cooling DP and the
    repump open; it is cut the same way, so each atom is recooled for exactly
    the recooling time.

    The duration is the shared t_recooling (t_recooling_two_nodes), as for
    the pumping stages in record_CW_optical_pumping_both_nodes.

    Leaves the coils at the PGC setting, like the original. The caller sets
    the OP coils again afterwards (Two_nodes_HOM_experiment does, as
    single_photon_experiment_9 does after record_recooling).
    """
    t_recooling = self.t_recooling

    with self.core_dma.record("recooling_both_nodes"):
        ### PGC coils on both nodes, loaded on one timestamp
        self.zotino0_Node1.write_dac(self.coil_channels_Node1[0], self.AZ_bottom_volts_PGC_Node1)
        self.zotino0_Node1.write_dac(self.coil_channels_Node1[1], -self.AZ_bottom_volts_PGC_Node1)
        self.zotino0_Node1.write_dac(self.coil_channels_Node1[2], self.AX_volts_PGC_Node1)
        self.zotino0_Node1.write_dac(self.coil_channels_Node1[3], self.AY_volts_PGC_Node1)
        self.zotino0_Node2.write_dac(self.coil_channels_Node2[0], self.AZ_bottom_volts_PGC_Node2)
        self.zotino0_Node2.write_dac(self.coil_channels_Node2[1], -self.AZ_bottom_volts_PGC_Node2)
        self.zotino0_Node2.write_dac(self.coil_channels_Node2[2], self.AX_volts_PGC_Node2)
        self.zotino0_Node2.write_dac(self.coil_channels_Node2[3], self.AY_volts_PGC_Node2)
        delay(1.5 * us)
        self.zotino0_Node1.load()
        self.zotino0_Node2.load()

        self.dds_cooling_DP_Node1.set(frequency=self.f_cooling_DP_PGC_Node1, amplitude=self.ampl_cooling_DP_PGC_Node1)
        self.dds_cooling_DP_Node2.set(frequency=self.f_cooling_DP_PGC_Node2, amplitude=self.ampl_cooling_DP_PGC_Node2)
        delay(0.4 * ms)  ### coil relaxation time

        ### fiber AOMs first: no light reaches an atom until its cooling DP opens
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
        delay(5 * us)

        ### RECOOLING STARTS HERE -- one timestamp, both nodes
        self.dds_cooling_DP_Node1.sw.on()
        self.dds_cooling_DP_Node2.sw.on()
        self.ttl_repump_switch_Node1.off()
        self.ttl_repump_switch_Node2.off()

        delay(t_recooling)

        ### ... AND ENDS HERE -- one timestamp, both nodes
        self.dds_cooling_DP_Node1.sw.off()
        self.dds_cooling_DP_Node2.sw.off()
        self.ttl_repump_switch_Node1.on()
        self.ttl_repump_switch_Node2.on()
        delay(5 * us)

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
        if not self.PGC_and_RO_with_on_chip_beams_Node1:
            self.dds_AOM_A5_Node1.sw.off()
            self.dds_AOM_A6_Node1.sw.off()
        if not self.PGC_and_RO_with_on_chip_beams_Node2:
            self.dds_AOM_A5_Node2.sw.off()
            self.dds_AOM_A6_Node2.sw.off()


@kernel
def record_chopped_blow_away_both_nodes(self):
    """Record ONE DMA trace "chopped_blow_away_both_nodes": blow-away on both nodes.

    Two-node version of record_chopped_blow_away (experiment_functions.py).
    There the blow-away light is the cooling DP light through fiber AOM A6
    alone (A1-A5 off), held on while the FORT is chopped -- off for 35% of each
    t_BA_chop_period -- so the atom sees the light unshifted.

    The nodes disagree on nearly every number: t_blowaway (150 vs 30 us live),
    t_BA_chop_period (1.2 vs 1.0 us), the cooling frequency and the coils. So
    the two chop trains are NOT interleaved: Node1 is blown first, then Node2,
    each with exactly its single-node timing, in one trace whose timestamps
    only increase. That costs the second node's t_blowaway, and each atom just
    sits in its FORT, with its own light off, while the other node is blown.

    Not used by Two_nodes_HOM_experiment, which reads out no spin state; it is
    for the two-node state readout that follows the HOM measurement.
    """
    n_chop_Node1 = int(self.t_blowaway_Node1 / self.t_BA_chop_period_Node1 + 0.5)
    n_chop_Node2 = int(self.t_blowaway_Node2 / self.t_BA_chop_period_Node2 + 0.5)
    assert n_chop_Node1 >= 1, "t_blowaway_Node1 should be >= t_BA_chop_period_Node1"
    assert n_chop_Node2 >= 1, "t_blowaway_Node2 should be >= t_BA_chop_period_Node2"

    ### the FORT is off for 35% of each chop period, as in the single-node trace
    period_Node1_mu = self.core.seconds_to_mu(self.t_BA_chop_period_Node1)
    pulse_Node1_mu = self.core.seconds_to_mu(self.t_BA_chop_period_Node1 * 0.35)
    period_Node2_mu = self.core.seconds_to_mu(self.t_BA_chop_period_Node2)
    pulse_Node2_mu = self.core.seconds_to_mu(self.t_BA_chop_period_Node2 * 0.35)

    with self.core_dma.record("chopped_blow_away_both_nodes"):
        self.ttl_repump_switch_Node1.on()  ### Repump AOM OFF
        self.ttl_repump_switch_Node2.on()
        delay(5 * us)

        ### Blow-away coils on both nodes, loaded on one timestamp
        self.zotino0_Node1.write_dac(self.coil_channels_Node1[0], self.AZ_bottom_volts_blowaway_Node1)
        self.zotino0_Node1.write_dac(self.coil_channels_Node1[1], -self.AZ_bottom_volts_blowaway_Node1)
        self.zotino0_Node1.write_dac(self.coil_channels_Node1[2], self.AX_volts_blowaway_Node1)
        self.zotino0_Node1.write_dac(self.coil_channels_Node1[3], self.AY_volts_blowaway_Node1)
        self.zotino0_Node2.write_dac(self.coil_channels_Node2[0], self.AZ_bottom_volts_blowaway_Node2)
        self.zotino0_Node2.write_dac(self.coil_channels_Node2[1], -self.AZ_bottom_volts_blowaway_Node2)
        self.zotino0_Node2.write_dac(self.coil_channels_Node2[2], self.AX_volts_blowaway_Node2)
        self.zotino0_Node2.write_dac(self.coil_channels_Node2[3], self.AY_volts_blowaway_Node2)
        delay(1.5 * us)
        self.zotino0_Node1.load()
        self.zotino0_Node2.load()

        delay(1 * ms)  ### Physical coil-relaxation time

        ### Configure the blow-away light on both nodes, with everything but A6 off
        self.dds_cooling_DP_Node1.set(frequency=self.f_cooling_DP_blowaway_Node1, amplitude=self.ampl_cooling_DP_MOT_Node1)
        self.dds_cooling_DP_Node2.set(frequency=self.f_cooling_DP_blowaway_Node2, amplitude=self.ampl_cooling_DP_MOT_Node2)
        delay(5 * us)
        self.dds_AOM_A1_Node1.sw.off()
        self.dds_AOM_A2_Node1.sw.off()
        self.dds_AOM_A3_Node1.sw.off()
        self.dds_AOM_A4_Node1.sw.off()
        delay(5 * us)
        self.dds_AOM_A1_Node2.sw.off()
        self.dds_AOM_A2_Node2.sw.off()
        self.dds_AOM_A3_Node2.sw.off()
        self.dds_AOM_A4_Node2.sw.off()
        delay(5 * us)
        self.dds_AOM_A5_Node1.sw.off()
        self.dds_AOM_A5_Node2.sw.off()
        ### Both nodes' blow-away light off until each node's own window,
        ### whatever the caller left on: the waiting node must not be blown
        ### (unchopped) while the other one is.
        self.dds_cooling_DP_Node1.sw.off()
        self.dds_cooling_DP_Node2.sw.off()
        delay(5 * us)
        self.dds_AOM_A6_Node1.sw.off()
        self.dds_AOM_A6_Node2.sw.off()
        delay(5 * us)
        if not self.blowaway_light_off_Node1:
            self.dds_AOM_A6_Node1.set(frequency=self.AOM_A6_freq_Node1, amplitude=dB_to_V(-7.0))
        if not self.blowaway_light_off_Node2:
            self.dds_AOM_A6_Node2.set(frequency=self.AOM_A6_freq_Node2, amplitude=dB_to_V(-7.0))
        delay(5 * us)

        ### Node1: light on with the first FORT-off, chopped FORT, light off
        start_Node1 = now_mu()
        if not self.blowaway_light_off_Node1:
            self.dds_AOM_A6_Node1.sw.on()
            self.dds_cooling_DP_Node1.sw.on()
        for i in range(n_chop_Node1):
            at_mu(start_Node1 + i * period_Node1_mu)
            self.dds_FORT_Node1.sw.off()

            at_mu(start_Node1 + i * period_Node1_mu + pulse_Node1_mu)
            self.dds_FORT_Node1.sw.on()

        at_mu(start_Node1 + n_chop_Node1 * period_Node1_mu)
        self.dds_FORT_Node1.sw.on()
        self.dds_cooling_DP_Node1.sw.off()
        self.ttl_repump_switch_Node1.on()
        self.dds_AOM_A6_Node1.sw.off()
        delay(5 * us)

        ### Node2: the same, after Node1
        start_Node2 = now_mu()
        if not self.blowaway_light_off_Node2:
            self.dds_AOM_A6_Node2.sw.on()
            self.dds_cooling_DP_Node2.sw.on()
        for i in range(n_chop_Node2):
            at_mu(start_Node2 + i * period_Node2_mu)
            self.dds_FORT_Node2.sw.off()

            at_mu(start_Node2 + i * period_Node2_mu + pulse_Node2_mu)
            self.dds_FORT_Node2.sw.on()

        at_mu(start_Node2 + n_chop_Node2 * period_Node2_mu)
        self.dds_FORT_Node2.sw.on()
        self.dds_cooling_DP_Node2.sw.off()
        self.ttl_repump_switch_Node2.on()
        self.dds_AOM_A6_Node2.sw.off()
        delay(5 * us)

        ### Restore the normal RF settings, as the single-node trace does
        self.dds_cooling_DP_Node1.set(frequency=self.f_cooling_DP_RO_Node1, amplitude=self.ampl_cooling_DP_MOT_Node1)
        self.dds_cooling_DP_Node2.set(frequency=self.f_cooling_DP_RO_Node2, amplitude=self.ampl_cooling_DP_MOT_Node2)
        self.dds_AOM_A6_Node1.set(frequency=self.AOM_A6_freq_Node1, amplitude=self.stabilizer_AOM_A6_Node1.amplitude)
        self.dds_AOM_A6_Node2.set(frequency=self.AOM_A6_freq_Node2, amplitude=self.stabilizer_AOM_A6_Node2.amplitude)


@kernel
def record_excitation_and_photon_collection_both_nodes(self) -> TInt64:
    """Record ONE DMA trace "excitation_and_photon_collection_both_nodes".

    Two-node version of record_excitation_and_photon_collection
    (experiment_functions.py). In single-node mode, relative to t1 (5 us after
    the trace starts), each node's FORT goes off at t1 and back on at
    t1 + 50 + collection, its excitation pulse (ttl_GRIN2_switch) fires at
    t1 + t_excitation_offset_mu, and the SPCMs count from t1 +
    gate_start_offset_mu for t_photon_collection_time. Each node's numbers are
    calibrated in that frame (its FORT and excitation AOM latencies), so here
    each node plays exactly that block, Node2's shifted by

        D = t_Node2_excitation_delay_mu - t_Node2_rtio_offset_mu

    with two changes forced by sharing ONE gate:

    ONE GATE FOR BOTH PHOTONS. The four SPCMs are shared, so they are gated
    once, over the union of the two nodes' (shifted) gates:
        gate_start = min(G1, G2 + D),  gate_end = max(G1, G2 + D) + L
    where G = gate_start_offset_mu per node and L = t_photon_collection_time,
    the shared two-node setting (t_photon_collection_time_two_nodes).

    BOTH FORTS DARK FOR THE WHOLE GATE. Node k's FORT goes off at
    t1 + gate_start - Gk and back on at t1 + 50 + gate_end - Gk: its own
    single-node block, widened only as far as the shared gate overhangs its
    own. Each node keeps its calibrated lead from FORT-off to gate-open (Gk)
    and from FORT-on to gate-close (Gk - 50), so no FORT light from either node
    reaches the SPCMs while they count. The two nodes' FORT edges are G2 - G1
    apart in programmed time (140 mu live) at every D, which is what aligns
    them at the detectors. Both FORTs stay dark for 50 + gate length.

    THE SIGN AND SIZE OF D. Positive moves Node2 LATER. With both blocks
    started together, Node2's photons arrive ~143 ns after Node1's
    (single-node RIDs 39221 vs 38894/39511, 2026_Compare_Single_Photon_Arrival
    _Timing; this already includes the satellite's ~281 ns DRTIO counter
    offset), so coincidence needs D ~ -143. Confirm it with this trace and one
    node's excitation off (override p_GRIN2_excitation_NodeX = -150.0) by
    comparing the two arrival histograms. t_Node2_excitation_delay_mu is the
    delay itself, as its declaration asks of ported code: its stored 189 is the
    standalone handshake's zero point and must be changed before this runs.
    t_Node2_rtio_offset_mu (electrical, stored 0) is subtracted, as wherever
    else Node2 events are placed; if it is ever set to the measured ~281,
    t_Node2_excitation_delay_mu becomes the purely optical delay (~ +138).

    THE TRACE. Both ttl_exc0_switch open at the trace start, 5 us before t1,
    as in single-node mode -- not shifted by D, which could put Node2's before
    the trace start. Every other edge is placed with at_mu, at most four on one
    timestamp, and the recording ends on the last edge so a playback leaves
    the cursor there. The satellite's edges come in time order (one step back
    for D >= ~100); the master's step back 4-5 times per playback (the three
    repeated gate openings and the exc0 close), using at most 6 of its 8 SED
    lanes for |D| <= 400 ns (modelled, 2026-10-26: no sequence errors).

    Returns gate_end, in mu after t1: the caller reads the SPCMs with
    timestamp_mu(t1 + gate_end), where t1 = playback start + 5 us.
    """
    D = int64(self.t_Node2_excitation_delay_mu) - int64(self.t_Node2_rtio_offset_mu)
    t1 = self.core.seconds_to_mu(5 * us)

    ### each node's single-node numbers, Node2's excitation shifted by D
    excitation_Node1 = int64(self.t_excitation_offset_mu_Node1)
    excitation_Node2 = D + int64(self.t_excitation_offset_mu_Node2)
    pulse_Node1 = self.core.seconds_to_mu(self.t_excitation_pulse_Node1)
    pulse_Node2 = self.core.seconds_to_mu(self.t_excitation_pulse_Node2)
    gate_offset_Node1 = int64(self.gate_start_offset_mu_Node1)
    gate_offset_Node2 = int64(self.gate_start_offset_mu_Node2)
    ### one collection time for both nodes (t_photon_collection_time_two_nodes)
    collection_Node1 = self.core.seconds_to_mu(self.t_photon_collection_time)
    collection_Node2 = collection_Node1

    ### the shared gate: the union of the two (shifted) gates
    gate_start = gate_offset_Node1
    if gate_offset_Node2 + D < gate_start:
        gate_start = gate_offset_Node2 + D
    gate_end = gate_offset_Node1 + collection_Node1
    if gate_offset_Node2 + D + collection_Node2 > gate_end:
        gate_end = gate_offset_Node2 + D + collection_Node2

    ### each FORT dark around the shared gate, with its own calibrated leads
    FORT_off_Node1 = gate_start - gate_offset_Node1
    FORT_on_Node1 = 50 + gate_end - gate_offset_Node1
    FORT_off_Node2 = gate_start - gate_offset_Node2
    FORT_on_Node2 = 50 + gate_end - gate_offset_Node2

    ### excitation path closed 1 us after each pulse, as in single-node mode
    exc0_close_Node1 = excitation_Node1 + pulse_Node1 + 1000
    exc0_close_Node2 = excitation_Node2 + pulse_Node2 + 1000

    trace_end = gate_end
    if exc0_close_Node1 > trace_end:
        trace_end = exc0_close_Node1
    if exc0_close_Node2 > trace_end:
        trace_end = exc0_close_Node2
    if FORT_on_Node1 > trace_end:
        trace_end = FORT_on_Node1
    if FORT_on_Node2 > trace_end:
        trace_end = FORT_on_Node2

    ### Nothing may land before the exc0 openings at the trace start: that
    ### would be a rewind on both crates, and with a large negative D an event
    ### before the playback start. ~ -4.6 us of headroom.
    earliest = excitation_Node1
    if excitation_Node2 < earliest:
        earliest = excitation_Node2
    if FORT_off_Node1 < earliest:
        earliest = FORT_off_Node1
    if FORT_off_Node2 < earliest:
        earliest = FORT_off_Node2
    if t1 + earliest < 8:
        raise ValueError("Node2 excitation delay too negative for the 5 us lead in the excitation trace")

    with self.core_dma.record("excitation_and_photon_collection_both_nodes"):
        self.ttl_exc0_switch_Node1.off()
        self.ttl_exc0_switch_Node2.off()

        ### the excitation pulses, each in its own node's block
        at_mu(t1 + excitation_Node2)
        self.ttl_GRIN2_switch_Node2.off()  # turns on excitation
        at_mu(t1 + excitation_Node2 + pulse_Node2)
        self.ttl_GRIN2_switch_Node2.on()  # turns off excitation
        at_mu(t1 + excitation_Node1)
        self.ttl_GRIN2_switch_Node1.off()
        at_mu(t1 + excitation_Node1 + pulse_Node1)
        self.ttl_GRIN2_switch_Node1.on()

        ### both FORTs dark around the shared gate
        at_mu(t1 + FORT_off_Node2)
        self.dds_FORT_Node2.sw.off()
        at_mu(t1 + FORT_off_Node1)
        self.dds_FORT_Node1.sw.off()
        at_mu(t1 + FORT_on_Node2)
        self.dds_FORT_Node2.sw.on()
        at_mu(t1 + FORT_on_Node1)
        self.dds_FORT_Node1.sw.on()

        ### ONE gate on the four shared, master-local SPCMs
        at_mu(t1 + gate_start)
        self.ttl_SPCM0.gate_rising_mu(gate_end - gate_start)
        at_mu(t1 + gate_start)
        self.ttl_SPCM1.gate_rising_mu(gate_end - gate_start)
        at_mu(t1 + gate_start)
        self.ttl_SPCM0_OtherNode.gate_rising_mu(gate_end - gate_start)
        at_mu(t1 + gate_start)
        self.ttl_SPCM1_OtherNode.gate_rising_mu(gate_end - gate_start)

        ### close each excitation path
        at_mu(t1 + exc0_close_Node2)
        self.ttl_exc0_switch_Node2.on()
        at_mu(t1 + exc0_close_Node1)
        self.ttl_exc0_switch_Node1.on()

        at_mu(t1 + trace_end)

    return gate_end


@rpc
def set_780_waveplates_both_nodes(self):
    """Move both nodes' 780 HWP and QWP to target_780_HWP/QWP_NodeX.

    Host RPC that blocks until all four rotators have stopped (the K10CR1
    driver's move_to waits for the move), so the caller must resync the
    timeline afterwards. Base publishes one node-aware view of the shared
    k10cr1_ndsp controller per node in two-node mode, which turns '780_HWP'
    into '780_HWP_Node1' / '780_HWP_Node2'.
    """
    for node in ("Node1", "Node2"):
        controller = getattr(self, "k10cr1_ndsp_" + node)
        deg_to_pos = getattr(self, "deg_to_pos_" + node)
        for plate in ("780_HWP", "780_QWP"):
            target_deg = getattr(self, "target_" + plate + "_" + node)
            controller.move_to(int(target_deg * deg_to_pos), plate)


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
    ###
    ### WORTH IT ONLY BECAUSE THE BRANCHES CONTAIN delay(). The manual: "the
    ### statements in a parallel block are actually executed sequentially,
    ### only the RTIO events are SCHEDULED in parallel". Two switch calls with
    ### nothing between them already share a timestamp, so wrapping THOSE in
    ### parallel changes no timing and just spends CPU -- an at_mu, a now_mu,
    ### a compare and a select per statement.
    ###
    ### And the rewind is not free: wall clock advances across a branch even
    ### though the cursor is reset, so a branch that emits many events can put
    ### the NEXT branch into negative slack. Two or three events per branch,
    ### as here, is nowhere near that.
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
        max_rounds = self.max_loading_rounds
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

        ### Lower both FORTs to their science set points BEFORE PGC, as the
        ### single-node loaders do (load_MOT_and_FORT_until_atom_recycle).
        ### Until 2026-10-27 this loader ran PGC at the loading depth (5.8x
        ### science on Node1, 2.9x on Node2) and only first_shot stepped the
        ### FORT down. Only a fresh load needs it: a reused atom is already at
        ### science from the previous second shot. Feedback-corrected science
        ### amplitude (amplitudes[1]), never p_FORT_*.
        self.dds_FORT_Node1.set(frequency=self.f_FORT_Node1,
                                amplitude=self.stabilizer_FORT_Node1.amplitudes[1])
        self.dds_FORT_Node2.set(frequency=self.f_FORT_Node2,
                                amplitude=self.stabilizer_FORT_Node2.amplitudes[1])
        delay(5 * us)

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
def load_until_atom_in_both_nodes_alternating_shot_recycle(self):
    """
    Two-node combined loading, with a PER-NODE atom check.

    ###################################################################
    ### NOT WORKING -- DO NOT USE. NEEDS MORE MODIFICATION.         ###
    ###################################################################
    As of 2026-10-08 this loader does not load atoms. Use
    load_until_atom_in_both_nodes_recycle (the joint check) instead;
    that one works and is unchanged.

    What is known so far, from RID 38718 (this loader) against RID 38710
    (the joint loader, working):

    * The atom is never loaded -- it is NOT lost after loading. The
      readout saw 0 atoms in 100 shots on Node2 and 1 on Node1, while
      the loader reported 100/100 loaded. So the check passes on
      background and the sequence proceeds with empty traps.
    * Part of it is threshold calibration. A window with one node dark
      still collects the FULL common background (~39900 c/s of FORT
      scatter plus dark counts, not gated by the cooling DP this grid
      switches) and only HALF the MOT fluorescence (~16300 c/s total),
      so the per-node no-atom floor is ~48050 c/s -- which is also what
      you get by halving the joint 96100. Halving therefore puts the
      threshold exactly on the floor. See the derivation beside
      alternating_shot_threshold_for_loading_Node1 in
      ExperimentVariables_master_satellite_Node1.py.
    * Raising the thresholds was NOT sufficient: it still does not load.
      So there is at least one further fault beyond calibration.
    * Unverified suspicion, and the place to start next: the MOT beams
      (dds_AOM_A1..A6 per node), the FORTs and the coils may not be
      left in the right state by this loader's window grid. The grid
      toggles each node's cooling DP and repump every window, which the
      joint loader never does, and nothing here re-asserts the MOT AOMs
      or the coil DACs afterwards.
    * Separately, and true of BOTH two-node loaders: PGC runs at the
      LOADING FORT depth. All four single-node loaders drop the FORT to
      stabilizer_FORT.amplitudes[1] before PGC. Not the cause of the 0%
      above, but a real porting gap.

    Note for whoever picks this up: nothing calls this function any more, so
    it is NOT on a compiled path and the ARTIQ compiler no longer type-checks
    it. Re-point Two_nodes_alternating_shot_experiment at it before trusting
    that it still compiles, and put "Two_nodes_alternating_shot_experiment"
    back into TWO_NODE_ALTERNATING_EXPERIMENTS in applets_master_satellite.py
    so the per-node atom-check histograms come up again.

    Below is the original design intent, kept for whoever picks this up.
    ###################################################################

    Identical to load_until_atom_in_both_nodes_recycle except for the check.
    That one opens ONE gate of all four SPCMs and compares the sum against the
    joint two_atom_threshold_for_loading -- and because every detector sees
    both traps through the beamsplitter fan-out, that sum structurally cannot
    tell "one atom in each trap" from "one bright atom plus background".

    Here the check alternates instead: each window darkens the OTHER node, so
    the four shared counters see only the active node's trap, and each node
    gets its own number compared against its own threshold. Both must pass.

    UNDER MOT LIGHT, not readout light. The loader is already in its MOT state
    when the check runs and stays there; only which node's cooling DP + repump
    is open alternates. So the thresholds are the _for_loading pair (brighter),
    not the readout pair.

    COST, measured against RID 38710 (joint loader, median load 0.85 s ~ 84
    tries): a try goes from 10.1 ms to 2*n*(window + 3*pad) = 26 ms at the
    defaults, so expect the median to land near 2 s. On top of that each node's
    MOT light is off for half the grid, so loading runs at ~50% duty WHILE
    checking -- trapped atoms are not lost, since the FORT holds them
    throughout, but new loading pauses. Expect the upper end of 2-3x.

    If that is too slow, the cheapest lever is a joint pre-check before this
    one: 51000 + 45000 = 96000 > 89000, so a per-node pass ALWAYS implies a
    joint pass and the old 10 ms gate is a filter that cannot reject a real
    load.
    """

    ### Validate the window geometry BEFORE touching hardware, so a bad knob
    ### cannot leave the MOT on. All three are scannable from the dashboard
    ### and every event in the grid below is placed absolutely with at_mu, so
    ### a pad too small to cover the Node2 offset turns an at_mu into a
    ### BACKWARDS jump and underflows mid-check rather than merely running a
    ### short pad.
    if self.t_alternating_RO_pad <= 0.0:
        raise ValueError(
            "t_alternating_RO_pad must be > 0: it is the guard band that "
            "keeps one node's MOT light out of the other node's counted "
            "window, and the check places every event absolutely from it.")
    if self.t_alternating_RO_window <= 0.0:
        raise ValueError("t_alternating_RO_window must be > 0.")
    if (self.t_Node2_rtio_offset_mu
            >= self.core.seconds_to_mu(self.t_alternating_RO_pad)
            or -self.t_Node2_rtio_offset_mu
            >= self.core.seconds_to_mu(self.t_alternating_RO_pad)):
        raise ValueError(
            "|t_Node2_rtio_offset_mu| must be smaller than "
            "t_alternating_RO_pad, or the Node2 shift eats the guard band "
            "separating the two nodes' light.")

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
        max_rounds = self.max_loading_rounds

        ### The same window geometry the readout uses, so the two stages stay
        ### comparable and there is one set of knobs to tune.
        n_windows = self.n_alternating_RO_windows_per_node
        pad_mu = self.core.seconds_to_mu(self.t_alternating_RO_pad)
        window_mu = self.core.seconds_to_mu(self.t_alternating_RO_window)
        node2_offset_mu = self.t_Node2_rtio_offset_mu
        period_mu = window_mu + 3 * pad_mu

        ### Per-node counting time. NOT t_atom_check_time -- that is the joint
        ### gate's duration and does not apply here.
        exposure = n_windows * self.t_alternating_RO_window

        AllSPCMs_atom_check = 0
        AllSPCMs_atom_check_not_loaded = 0
        AllSPCMs_atom_check_Node1 = 0
        AllSPCMs_atom_check_Node2 = 0
        AllSPCMs_atom_check_not_loaded_Node1 = 0
        AllSPCMs_atom_check_not_loaded_Node2 = 0
        rounds = 0

        while not atom_loaded and rounds < max_rounds:
            rounds += 1
            ### try_n is reset EVERY round. In the legacy version it was set
            ### once outside the loop and reset only inside the laser-feedback
            ### branch, so with feedback off the inner loop stopped checking
            ### after max_tries and the outer loop spun doing nothing.
            try_n = 0

            while not atom_loaded and try_n < max_tries:
                ### ALTERNATING ATOM CHECK, under MOT light.
                ###
                ### Both nodes' MOT light is ON when we arrive here. Each
                ### window darkens the OTHER node so the four shared,
                ### master-local counters see only the active node's trap,
                ### then puts it back. Everything is placed absolutely from
                ### one t_grid, exactly as two_node_alternating_shot does.
                self.core.break_realtime()
                t_grid = now_mu()

                for i in range(2 * n_windows):
                    t_window = t_grid + i * period_mu
                    node1_window = (i % 2) == 0

                    ### dark the OTHER node for this window
                    if node1_window:
                        at_mu(t_window - node2_offset_mu)
                        self.ttl_repump_switch_Node2.on()
                        self.dds_cooling_DP_Node2.sw.off()
                    else:
                        at_mu(t_window)
                        self.ttl_repump_switch_Node1.on()
                        self.dds_cooling_DP_Node1.sw.off()

                    ### the counters are master-local, so never offset
                    at_mu(t_window + pad_mu)
                    with parallel:
                        self.ttl_SPCM0_counter.gate_rising_mu(window_mu)
                        self.ttl_SPCM1_counter.gate_rising_mu(window_mu)
                        self.ttl_SPCM0_OtherNode_counter.gate_rising_mu(window_mu)
                        self.ttl_SPCM1_OtherNode_counter.gate_rising_mu(window_mu)

                    ### TRAILING PAD, and it is not cosmetic. gate_rising_mu
                    ### emits a CLOSE event as well as an open one, so all
                    ### four counters close at t_window + pad + window.
                    ### Restoring the light at that same instant would put SIX
                    ### events on one timestamp -- past the four-event rule,
                    ### and two short of the eight SED FIFOs beyond which
                    ### events are silently DISCARDED.
                    if node1_window:
                        at_mu(t_window + 2 * pad_mu + window_mu - node2_offset_mu)
                        self.ttl_repump_switch_Node2.off()
                        self.dds_cooling_DP_Node2.sw.on()
                    else:
                        at_mu(t_window + 2 * pad_mu + window_mu)
                        self.ttl_repump_switch_Node1.off()
                        self.dds_cooling_DP_Node1.sw.on()

                at_mu(t_grid + 2 * n_windows * period_mu)
                delay(0.1 * ms)

                ### Deferred fetches, in emission order. The parity of the
                ### FETCH index reproduces the parity of the gate, which is
                ### what attributes the counts to a node. 2*n totals per
                ### counter against the 64-deep input FIFO, so this is safe to
                ### n = 32; past that the failure is a loud RTIOOverflow.
                AllSPCMs_atom_check_Node1 = 0
                AllSPCMs_atom_check_Node2 = 0
                for i in range(2 * n_windows):
                    window_count = (
                        self.ttl_SPCM0_counter.fetch_count()
                        + self.ttl_SPCM1_counter.fetch_count()
                        + self.ttl_SPCM0_OtherNode_counter.fetch_count()
                        + self.ttl_SPCM1_OtherNode_counter.fetch_count())
                    if (i % 2) == 0:
                        AllSPCMs_atom_check_Node1 += window_count
                    else:
                        AllSPCMs_atom_check_Node2 += window_count

                ### kept so the joint histogram stays comparable with every
                ### run taken before this function existed
                AllSPCMs_atom_check = (AllSPCMs_atom_check_Node1
                                       + AllSPCMs_atom_check_Node2)

                try_n += 1
                if try_n == 1:
                    AllSPCMs_atom_check_not_loaded = AllSPCMs_atom_check
                    AllSPCMs_atom_check_not_loaded_Node1 = AllSPCMs_atom_check_Node1
                    AllSPCMs_atom_check_not_loaded_Node2 = AllSPCMs_atom_check_Node2

                ### BOTH nodes, each against its own threshold. This is the
                ### whole point: the joint sum could be satisfied by one
                ### bright trap, this cannot.
                ###
                ### The thresholds are NOT the joint one halved, and getting
                ### that wrong is silent -- it makes the loader report 100%
                ### loading with 0% retention rather than raising. A window
                ### with one node dark still collects the FULL common
                ### background (FORT scatter + dark counts, ~39900 c/s, not
                ### gated by the cooling DP this grid switches) and only HALF
                ### the MOT fluorescence (~16300 c/s total), so the per-node
                ### no-atom floor is ~48050 c/s while half the joint value is
                ### ~48050 too -- i.e. halving puts the threshold exactly on
                ### the floor. An atom adds ~25100 c/s; the valley is ~60600.
                ### Full derivation beside the Variable in
                ### ExperimentVariables_master_satellite_Node1.py.
                ###
                ### For the same reason two_atom_threshold_for_loading (89000)
                ### must NOT be reused on AllSPCMs_atom_check below: the sum of
                ### two single-lit windows double-counts the background, and
                ### 91% of RID 38718's sums cleared 89000 with no atom present.
                if (AllSPCMs_atom_check_Node1/exposure
                        > self.alternating_shot_threshold_for_loading_Node1
                        and AllSPCMs_atom_check_Node2/exposure
                        > self.alternating_shot_threshold_for_loading_Node2):
                    atom_loaded = True
                    t_after_atom = now_mu()

            if atom_loaded:
                self.set_dataset("time_without_atom", 0.0, broadcast=True)
                ### Two values per loaded round, as the joint loader does:
                ### the count that crossed, and the first try's count that did
                ### not. That pairing is what makes the histogram bimodal and
                ### lets the threshold be read off the valley.
                self.append_to_dataset("AllSPCMs_atom_check_in_loading", AllSPCMs_atom_check)
                self.append_to_dataset("AllSPCMs_atom_check_in_loading", AllSPCMs_atom_check_not_loaded)
                self.append_to_dataset("AllSPCMs_alternating_RO_atom_check_in_loading_Node1", AllSPCMs_atom_check_Node1)
                self.append_to_dataset("AllSPCMs_alternating_RO_atom_check_in_loading_Node1", AllSPCMs_atom_check_not_loaded_Node1)
                self.append_to_dataset("AllSPCMs_alternating_RO_atom_check_in_loading_Node2", AllSPCMs_atom_check_Node2)
                self.append_to_dataset("AllSPCMs_alternating_RO_atom_check_in_loading_Node2", AllSPCMs_atom_check_not_loaded_Node2)
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
                "alternating_shot_threshold_for_loading_Node1/_Node2 against "
                "the per-node atom-check histograms.")

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

    ### The pads are the only headroom this loop has, and every knob that
    ### eats them is scannable from the dashboard. Every event below is placed
    ### ABSOLUTELY with at_mu, so a pad too small to cover the Node2 offset
    ### turns an at_mu into a BACKWARDS jump and the kernel underflows
    ### mid-loop rather than merely running a short pad; a pad of 0 collapses
    ### the window structure outright. Fail here, naming the knob.
    if self.t_alternating_RO_pad <= 0.0:
        raise ValueError(
            "t_alternating_RO_pad must be > 0: it is the guard band that "
            "keeps one node's readout light out of the other node's counted "
            "window, and the loop places every event absolutely from it."
        )
    if self.t_alternating_RO_window <= 0.0:
        raise ValueError("t_alternating_RO_window must be > 0.")

    n_windows = self.n_alternating_RO_windows_per_node
    pad_mu = self.core.seconds_to_mu(self.t_alternating_RO_pad)
    window_mu = self.core.seconds_to_mu(self.t_alternating_RO_window)
    node2_offset_mu = self.t_Node2_rtio_offset_mu
    period_mu = window_mu + 3 * pad_mu

    ### One pad now sits BEFORE the gate and one AFTER it, leaving the third
    ### as the dark band between nodes. Node2's events shift by the RTIO
    ### offset, eating that band from one side, so the offset must stay
    ### smaller than a pad or Node2's light can still be on when Node1's
    ### window opens.
    if node2_offset_mu >= pad_mu or -node2_offset_mu >= pad_mu:
        raise ValueError(
            "|t_Node2_rtio_offset_mu| must be smaller than "
            "t_alternating_RO_pad, or the Node2 shift eats the guard band "
            "separating the two nodes' readout light."
        )

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

        ### TRAILING PAD, and it is not cosmetic. gate_rising_mu emits a
        ### CLOSE event as well as an open one, so all four counters close at
        ### t_window + pad + window. Switching the light off at that same
        ### instant put SIX events on one timestamp -- past the four-event
        ### rule, and only two short of the eight SED FIFOs beyond which
        ### events are silently DISCARDED. One pad later splits that into
        ### four and two, leaves the 3-pad period untouched, and guarantees
        ### the whole counted window is lit instead of racing the light-off.
        if alice_window:
            at_mu(t_window + 2 * pad_mu + window_mu)
            self.dds_cooling_DP_Node1.sw.off()
            self.ttl_repump_switch_Node1.on()
        else:
            at_mu(t_window + 2 * pad_mu + window_mu - node2_offset_mu)
            self.dds_cooling_DP_Node2.sw.off()
            self.ttl_repump_switch_Node2.on()

    at_mu(t0 + 2 * n_windows * period_mu)
    delay(1 * ms)

    ### deferred fetches, in emission order
    self.AllSPCMs_alternating_RO_Node1 = 0
    self.AllSPCMs_alternating_RO_Node2 = 0
    for i in range(2 * n_windows):
        window_count = (self.ttl_SPCM0_counter.fetch_count()
                        + self.ttl_SPCM1_counter.fetch_count()
                        + self.ttl_SPCM0_OtherNode_counter.fetch_count()
                        + self.ttl_SPCM1_OtherNode_counter.fetch_count())
        if (i % 2) == 0:
            self.AllSPCMs_alternating_RO_Node1 += window_count
        else:
            self.AllSPCMs_alternating_RO_Node2 += window_count
    self.core.break_realtime()


@kernel
def two_node_alternating_windows(self):
    """Readout optics plus one alternating window grid. Returns (Node1, Node2).

    The shared primitive behind alternating_first_shot and
    alternating_second_shot. Factored rather than copied for two reasons, and
    the second one is the important one:

      * This grid has already needed two safety fixes -- the trailing pad for
        gate_rising_mu's CLOSE event, and the |offset| < pad validation -- and
        there are already two inline copies of it in this file. Copies three
        and four would diverge from those fixes silently.
      * RETENTION IS A RATIO, so the two shots must see IDENTICAL optics or the
        ratio means nothing. That is not hypothetical here: first_shot once
        lacked the A5/A6 else-branches that second_shot had, so the first
        readout ran on six beams and the second on four (see the long comment
        in first_shot). One function for both shots makes that class of
        asymmetry impossible rather than merely fixed.

    Deliberately does NOT write any dataset or attribute. The caller names what
    it measured, which is what lets one grid serve both shots.

    LIMIT ON n_alternating_RO_windows_per_node, unchanged from
    two_node_alternating_shot: the fetches are deferred, so each counter queues
    2 * n totals. At the default n = 10 that is 20 per counter against the
    64-deep Kasli input FIFO; past n = 32 it overflows. The failure is loud
    (RTIOOverflow, not a wrong number).
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

    ### Every event below is placed ABSOLUTELY with at_mu, so a pad too small
    ### to cover the Node2 offset turns an at_mu into a BACKWARDS jump and the
    ### kernel underflows mid-loop; a pad of 0 collapses the window structure
    ### outright. Both knobs are scannable from the dashboard, so fail here,
    ### naming the knob.
    if self.t_alternating_RO_pad <= 0.0:
        raise ValueError(
            "t_alternating_RO_pad must be > 0: it is the guard band that "
            "keeps one node's readout light out of the other node's counted "
            "window, and the loop places every event absolutely from it."
        )
    if self.t_alternating_RO_window <= 0.0:
        raise ValueError("t_alternating_RO_window must be > 0.")

    n_windows = self.n_alternating_RO_windows_per_node
    pad_mu = self.core.seconds_to_mu(self.t_alternating_RO_pad)
    window_mu = self.core.seconds_to_mu(self.t_alternating_RO_window)
    node2_offset_mu = self.t_Node2_rtio_offset_mu
    period_mu = window_mu + 3 * pad_mu

    if node2_offset_mu >= pad_mu or -node2_offset_mu >= pad_mu:
        raise ValueError(
            "|t_Node2_rtio_offset_mu| must be smaller than "
            "t_alternating_RO_pad, or the Node2 shift eats the guard band "
            "separating the two nodes' readout light."
        )

    self.core.break_realtime()
    t0 = now_mu()

    for i in range(2 * n_windows):
        t_window = t0 + i * period_mu
        node1_window = (i % 2) == 0

        ### only the active node turns on readout light
        if node1_window:
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

        ### TRAILING PAD, and it is not cosmetic. gate_rising_mu emits a CLOSE
        ### event as well as an open one, so all four counters close at
        ### t_window + pad + window. Switching the light off at that same
        ### instant puts SIX events on one timestamp -- past the four-event
        ### rule, and two short of the eight SED FIFOs beyond which events are
        ### silently DISCARDED.
        if node1_window:
            at_mu(t_window + 2 * pad_mu + window_mu)
            self.dds_cooling_DP_Node1.sw.off()
            self.ttl_repump_switch_Node1.on()
        else:
            at_mu(t_window + 2 * pad_mu + window_mu - node2_offset_mu)
            self.dds_cooling_DP_Node2.sw.off()
            self.ttl_repump_switch_Node2.on()

    at_mu(t0 + 2 * n_windows * period_mu)
    delay(1 * ms)

    ### deferred fetches, in emission order
    node1_counts = 0
    node2_counts = 0
    for i in range(2 * n_windows):
        window_count = (self.ttl_SPCM0_counter.fetch_count()
                        + self.ttl_SPCM1_counter.fetch_count()
                        + self.ttl_SPCM0_OtherNode_counter.fetch_count()
                        + self.ttl_SPCM1_OtherNode_counter.fetch_count())
        if (i % 2) == 0:
            node1_counts += window_count
        else:
            node2_counts += window_count
    self.core.break_realtime()

    return node1_counts, node2_counts


@kernel
def alternating_first_shot(self):
    """Node-resolved FIRST shot: one alternating grid, two numbers.

    The alternating counterpart of first_shot. Where that opens one gate of
    all four SPCMs and produces a single joint AllSPCMs_RO1 -- which cannot
    tell one atom in each trap from one bright atom plus background, because
    every detector sees both traps through the beamsplitter fan-out -- this
    separates the two nodes in time and produces one number per node.

    Writes only the attributes. The experiment appends them once both shots
    have run, so a measurement never lands half-recorded.
    """
    node1_counts, node2_counts = two_node_alternating_windows(self)
    self.AllSPCMs_alternating_RO1_Node1 = node1_counts
    self.AllSPCMs_alternating_RO1_Node2 = node2_counts


@kernel
def alternating_second_shot(self):
    """Node-resolved SECOND shot, for per-node retention.

    Identical optics and identical grid to alternating_first_shot -- the same
    function provides both -- so RO2/RO1 per node is a survival probability
    and not an artefact of the two shots being measured differently.

    Compare each node against alternating_shot_threshold_NodeX (counts/s),
    which was declared for exactly this and is read by nothing else yet.
    """
    node1_counts, node2_counts = two_node_alternating_windows(self)
    self.AllSPCMs_alternating_RO2_Node1 = node1_counts
    self.AllSPCMs_alternating_RO2_Node2 = node2_counts


@kernel
def end_measurement(self):
    """
    End the measurement by setting datasets and deciding whether to increment
    the measurement index.

    Copied from experiment_functions.end_measurement, with ONE difference: the
    magnetometer block is gone. It reads

        if self.monitor_magnetometer_in_end_measurement:
            if self.which_node == "alice": measure_Magnetometer_X_Y_Z(self)
            else:                          measure_Magnetometer_X_Y_Z(self)

    and although a two-node run must not take it (monitor_magnetometer_in_end
    _measurement is one shared global since 2026-10-27, and a two-node
    sequence samples no magnetometers -- which is why the block is gone), the
    ARTIQ
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
    self.append_to_dataset('AllSPCMs_alternating_RO_Node1_current_iteration', self.AllSPCMs_alternating_RO_Node1)
    self.append_to_dataset('AllSPCMs_alternating_RO_Node2_current_iteration', self.AllSPCMs_alternating_RO_Node2)
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

    self.append_to_dataset('AllSPCMs_alternating_RO_Node1', self.AllSPCMs_alternating_RO_Node1)
    self.append_to_dataset('AllSPCMs_alternating_RO_Node2', self.AllSPCMs_alternating_RO_Node2)
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

    TWO ADVANCE CRITERIA, selected by the global advance_if_atom_in_first_shot.
    This replaced a separate Two_nodes_atom_loading_2_experiment on 2026-10-12:
    the two differed only in this block, and one flag is cheaper to keep in
    step than two copies of a 100-line sequence.

      False (default, historical behaviour)
        Advance as soon as the LOADING CHECK says both traps are occupied.
        This is the only setting that MEASURES loading probability, because
        n_atom_loaded_per_iteration then counts one load per measurement.

      True
        Advance only once the first shot really saw two atoms, retrying
        load+first_shot up to max_loading_rounds times and recording
        the attempt count in n_first_shot_confirm_tries (and, since
        2026-10-27, the time from the first attempt to the passing shot in
        Atom_loading_time_first_shot_alive). Gives a FIXED number
        of confirmed two-atom events per scan point, which matters when the
        loading rate itself depends on the scanned parameter and would
        otherwise widen the error bars unevenly across the scan.

        Note what this does NOT fix: plot_retention_and_loading.py already
        post-selects, computing retention among shots with RO1 above
        threshold. So "only count when there was really an atom" is already
        how retention is defined -- the flag buys constant statistics, not
        correctness.

        On exhaustion the measurement is KEPT and the sequence carries on
        through the second shot and end_measurement, so a scan point always
        completes. Whether a given measurement was confirmed or merely
        exhausted is recoverable from the recorded AllSPCMs_RO1 itself.

    Both branches compile either way, which is the point of the flag living in
    a declared variable rather than arriving through an override dictionary: a
    missing override fails to COMPILE, not to run.
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

        if self.advance_if_atom_in_first_shot:
            ### Keep a measurement only once the FIRST SHOT really saw two
            ### atoms. Safe to retry because first_shot writes NO datasets --
            ### every append happens in end_measurement -- so a rejected
            ### attempt leaves no trace in any histogram or h5.
            ### Start of this measurement's whole load-until-alive span, for
            ### Atom_loading_time_first_shot_alive below. An elapsed-time
            ### stamp only, never an at_mu anchor.
            t_first_shot_alive_start = now_mu()
            atom_confirmed = False
            confirm_tries = 0
            while (not atom_confirmed
                   and confirm_tries < self.max_loading_rounds):
                confirm_tries += 1

                ### The loader's own already-loaded pre-check reads this. On a
                ### retry it is stale -- the shot below just said the trap is
                ### empty -- and leaving it set makes the loader return
                ### without loading anything, so this loop would spin without
                ### ever attempting to load again.
                self.AllSPCMs_RO2 = 0

                load_until_atom_in_both_nodes_recycle(self)
                delay(1 * ms)
                first_shot(self)
                delay(1 * ms)

                ### counts/s on BOTH sides. AllSPCMs_RO1 is raw counts over
                ### t_SPCM_first_shot and two_atom_threshold is a rate, so
                ### comparing them directly is ~390 > 39000 and never true.
                if (self.AllSPCMs_RO1 / self.t_SPCM_first_shot
                        > self.two_atom_threshold):
                    atom_confirmed = True
                else:
                    ### the rejected attempt's host RPCs advanced the wall
                    ### clock while the cursor stood still
                    self.core.break_realtime()

            ### The time it really took to get this measurement going: every
            ### load attempt, its PGC and first shot, until one shot saw both
            ### atoms (or max_loading_rounds ran out). Atom_loading_time has
            ### each load's MOT time only, so with the confirm on it leaves out
            ### the rejected loads and everything between a load and its shot.
            self.append_to_dataset(
                'Atom_loading_time_first_shot_alive',
                self.core.mu_to_seconds(now_mu() - t_first_shot_alive_start))

            ### Loading probability is 1 / mean(this series): the confirm loop
            ### makes n_atom_loaded_per_iteration count loads across attempts
            ### rather than measurements, so this is where that information
            ### goes. Not appended in the else branch below, so the series
            ### exists only for runs that actually used the flag.
            self.append_to_dataset('n_first_shot_confirm_tries', confirm_tries)
            self.core.break_realtime()

        else:
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
def Two_nodes_alternating_retention_experiment(self):
    """Per-node retention: BOTH shots node-resolved.

    Two_nodes_alternating_shot_experiment keeps the joint first and second
    shots and adds one alternating shot between them, so the retention it
    reports is still the joint RO2/RO1. This replaces both shots with
    alternating ones, which makes retention a per-node number:

        retention(NodeX) = fraction with RO2_NodeX above threshold
                           among those with RO1_NodeX above threshold

    compared against alternating_shot_threshold_Node1 / _Node2 (20000 / 21000
    counts/s), which were declared for exactly this and are read by nothing
    else yet.

    WHAT THIS DOES NOT PRODUCE, and it matters when you open the dashboard:
    the joint first_shot and second_shot never run here, so AllSPCMs_RO1 and
    AllSPCMs_RO2 stay at their seeded ZERO for the whole run. The existing
    "retention and loading AllSPCMs" applet will therefore read 0% -- it is not
    broken, it is being shown a dataset this experiment does not write. Read
    the four per-node histograms instead.

    Faking the joint numbers from the per-node sum was considered and
    rejected: a window with one node dark still carries the FULL common
    background, so the sum of two single-lit windows is the joint value PLUS
    one extra copy of the background and would read high against any joint
    threshold. That arithmetic is what made the per-node atom check fail --
    see the derivation beside alternating_shot_threshold_for_loading_Node1 in
    ExperimentVariables_master_satellite_Node1.py.

    The joint path is untouched: first_shot, second_shot,
    two_node_alternating_shot and both existing two-node experiments are
    unchanged.
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

        ### the JOINT loader, which works. The per-node loading check is a
        ### separate and currently broken path -- see
        ### load_until_atom_in_both_nodes_alternating_shot_recycle.
        load_until_atom_in_both_nodes_recycle(self)

        delay(1 * ms)

        alternating_first_shot(self)
        delay(1 * ms)

        if self.t_recooling_after_first_shot_Node1 > 0 or self.t_recooling_after_first_shot_Node2 > 0:
            recooling_after_first_shot(self)

        ### the readout leaves the fiber AOMs on, exactly as first_shot does
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

        ### t_FORT_drop is per node and is the knob retention is measured
        ### against, so each node gets exactly its own; `with parallel` resyncs
        ### the cursor past the longer drop.
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

        alternating_second_shot(self)

        ### Appended HERE, not inside the shots and not in end_measurement.
        ### Not in the shots: a measurement must never land half-recorded, and
        ### RO1 appended with RO2 missing would silently misalign every
        ### downstream zip of the two series. Not in end_measurement: that is
        ### shared with the joint experiments, and appending there would add
        ### four columns of zeros to every joint run's h5.
        self.append_to_dataset('AllSPCMs_alternating_RO1_Node1',
                               self.AllSPCMs_alternating_RO1_Node1)
        self.append_to_dataset('AllSPCMs_alternating_RO1_Node2',
                               self.AllSPCMs_alternating_RO1_Node2)
        self.append_to_dataset('AllSPCMs_alternating_RO2_Node1',
                               self.AllSPCMs_alternating_RO2_Node1)
        self.append_to_dataset('AllSPCMs_alternating_RO2_Node2',
                               self.AllSPCMs_alternating_RO2_Node2)
        self.append_to_dataset('AllSPCMs_alternating_RO1_Node1_current_iteration',
                               self.AllSPCMs_alternating_RO1_Node1)
        self.append_to_dataset('AllSPCMs_alternating_RO1_Node2_current_iteration',
                               self.AllSPCMs_alternating_RO1_Node2)
        self.append_to_dataset('AllSPCMs_alternating_RO2_Node1_current_iteration',
                               self.AllSPCMs_alternating_RO2_Node1)
        self.append_to_dataset('AllSPCMs_alternating_RO2_Node2_current_iteration',
                               self.AllSPCMs_alternating_RO2_Node2)
        self.core.break_realtime()

        end_measurement(self)

    self.append_to_dataset('n_feedback_per_iteration', self.n_feedback_per_iteration)
    self.append_to_dataset('n_atom_loaded_per_iteration', self.n_atom_loaded_per_iteration)


@kernel
def Two_nodes_HOM_experiment(self):
    """Hong-Ou-Mandel interference between single photons from Node1 and Node2.

    Master-satellite port of the standalone Two_node_single_photon_2_experiment
    (experiment_functions_obsolete.py), built on single_photon_experiment_9's
    loop and data format. Both nodes run ONE loop on one timeline, so the
    standalone TTL handshake is gone.

    1. Both nodes' 780 waveplates to target_780_HWP/QWP_NodeX (once per scan
       point, before anything else)
    2. Feedback, then the DMA traces are recorded
    3. Per measurement:
       a. load an atom in both nodes; first shot
       b. excitation loop, both nodes together, until either atom is lost:
            optical pumping           record_CW_optical_pumping_both_nodes
            n_excitation_attempts x   record_excitation_and_photon_collection
                                      _both_nodes, first click of each SPCM
                                      timestamped
            recooling every recool_every_n_OP cycles  record_recooling_both_nodes
            atom check every atom_check_every_n cycles: both atoms'
              fluorescence on the four SPCMs against two_atom_threshold
       c. second shot, end_measurement, the photon rows saved

    DATA. The compact format of single_photon_experiment_9, so the existing
    analysis code applies: one row per attempt in which any SPCM clicked --
    SPCM0/SPCM1/SPCM0_OtherNode/SPCM1_OtherNode_SinglePhoton_reduced_tStamps
    (first click, seconds, -1 if none) and reference_tStamps_t1 (that
    attempt's t1); per measurement n_photon_events (rows) and
    n_excitation_cycles (loop index at exit); AllSPCMs_RO_atom_check for every
    passed check. An HOM coincidence is a row with clicks on two SPCMs. The
    arrival time is click - t1, with the same t1 as in the single-node
    experiments.

    THE NODE2 DELAY is t_Node2_excitation_delay_mu, used directly as the delay
    (see record_excitation_and_photon_collection_both_nodes): about -143 to
    start with. Its stored 189 is the standalone zero point. Scanning it is
    the HOM scan.

    ONE NODE ONLY. To see each node's photon arrival in exactly this sequence,
    turn the other node's excitation off by override, p_GRIN2_excitation_Node2
    = -150.0 (or _Node1): that GRIN2 DDS then gets no RF, while everything
    else, the dark node's FORT switching included, runs unchanged.

    SHARED SETTINGS. The loop shape (n_excitation_attempts,
    max_excitation_cycles, atom_check_every_n, recool_every_n_OP) and the
    shared durations (t_pumping, t_depumping, t_OP_coil_relax_time,
    t_recooling, t_SPCM_recool_and_shot, t_photon_collection_time) are the
    two-node settings <x>_two_nodes, read here under the bare names. The
    per-node values single-node experiments use are not read. Override one
    for a run with the bare name in ONE node's dictionary, e.g.
    {'t_photon_collection_time': 150*ns}.

    The OP coils are set right after each recooling and each passed atom
    check, with 0.4 ms to settle, before the pumping trace sets them again --
    as in single_photon_experiment_9.
    """

    self.core.reset()
    delay(1 * ms)

    ### The shared loop shape (n_excitation_attempts_two_nodes etc.)
    n_excitation_attempts = self.n_excitation_attempts
    max_excitation_cycles = self.max_excitation_cycles
    atom_check_every_n = self.atom_check_every_n
    recool_every_n_OP = self.recool_every_n_OP

    ### One readout gate for both atoms in the atom check
    t_atom_check = self.t_SPCM_recool_and_shot
    do_pumping = self.t_pumping > 0.0
    do_recooling = self.t_recooling > 0.0

    ### photon rows kept per measurement, as in single_photon_experiment_9
    short_timestamps_length = int(max_excitation_cycles * 0.5)

    ### Both nodes' 780 waveplates. A host RPC that returns only once the
    ### rotators have stopped, so the timeline is resynced after it.
    set_780_waveplates_both_nodes(self)
    self.core.break_realtime()
    delay(10 * ms)

    self.n_feedback_per_iteration = 2  ### see Two_nodes_atom_loading_experiment
    self.n_atom_loaded_per_iteration = 0

    if self.enable_laser_feedback:
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP_Node1.set(frequency=self.f_cooling_DP_MOT_Node1,
                                      amplitude=self.ampl_cooling_DP_MOT_Node1)
        self.dds_cooling_DP_Node2.set(frequency=self.f_cooling_DP_MOT_Node2,
                                      amplitude=self.ampl_cooling_DP_MOT_Node2)
        run_feedback_and_record_FORT_MM_power(self)

    ### Recorded AFTER the feedback: the traces bake in the A5/A6 amplitudes
    ### they restore, and the excitation trace bakes in the Node2 delay and the
    ### shared gate of this scan point.
    record_CW_optical_pumping_both_nodes(self)
    record_recooling_both_nodes(self)
    gate_end_mu = record_excitation_and_photon_collection_both_nodes(self)
    ### In the log, so a run with the wrong delay (e.g. the stored 189) shows it
    self.print_async("Two_nodes_HOM_experiment: Node2 delay D =",
                     int64(self.t_Node2_excitation_delay_mu) - int64(self.t_Node2_rtio_offset_mu),
                     "mu; shared photon gate closes at t1 +", gate_end_mu, "mu")
    OP_handle = self.core_dma.get_handle("CW_optical_pumping_both_nodes")
    recooling_handle = self.core_dma.get_handle("recooling_both_nodes")
    excitation_handle = self.core_dma.get_handle("excitation_and_photon_collection_both_nodes")
    self.core.break_realtime()

    self.measurement = 0
    while self.measurement < self.n_measurements:

        AllSPCMs_RO_atom_check_array = [0] * int(max_excitation_cycles / atom_check_every_n)
        tStamps_t1 = [0.0] * short_timestamps_length
        SPCM0_timestamps = [-1.0] * short_timestamps_length
        SPCM1_timestamps = [-1.0] * short_timestamps_length
        SPCM0_OtherNode_timestamps = [-1.0] * short_timestamps_length
        SPCM1_OtherNode_timestamps = [-1.0] * short_timestamps_length
        Any_SPCM_click_counter = 0

        self.core.break_realtime()
        self.ttl_exc0_switch_Node1.on()  # turns off the excitation
        self.ttl_exc0_switch_Node2.on()
        delay(10 * ms)

        if self.advance_if_atom_in_first_shot:
            ### Keep a measurement only once the first shot really saw two
            ### atoms, as in Two_nodes_atom_loading_experiment.
            ### Start of this measurement's whole load-until-alive span, for
            ### Atom_loading_time_first_shot_alive below. An elapsed-time
            ### stamp only, never an at_mu anchor.
            t_first_shot_alive_start = now_mu()
            atom_confirmed = False
            confirm_tries = 0
            while (not atom_confirmed
                   and confirm_tries < self.max_loading_rounds):
                confirm_tries += 1

                load_until_atom_in_both_nodes_recycle(self)
                delay(1 * ms)
                first_shot(self)
                delay(1 * ms)

                if (self.AllSPCMs_RO1 / self.t_SPCM_first_shot
                        > self.two_atom_threshold):
                    atom_confirmed = True
                else:
                    ### The loader skips loading while the last second shot
                    ### saw both atoms. The shot above just said otherwise,
                    ### so clear it -- only here, as single_photon_experiment_9
                    ### does: clearing it before the first try would defeat
                    ### the atom reuse and force a full load every measurement.
                    self.AllSPCMs_RO2 = 0
                    self.core.break_realtime()

            ### The time it really took to get this measurement going: every
            ### load attempt, its PGC and first shot, until one shot saw both
            ### atoms (or max_loading_rounds ran out). Atom_loading_time has
            ### each load's MOT time only, so with the confirm on it leaves out
            ### the rejected loads and everything between a load and its shot.
            self.append_to_dataset(
                'Atom_loading_time_first_shot_alive',
                self.core.mu_to_seconds(now_mu() - t_first_shot_alive_start))
            self.append_to_dataset('n_first_shot_confirm_tries', confirm_tries)
            self.core.break_realtime()

        else:
            load_until_atom_in_both_nodes_recycle(self)

            delay(1 * ms)

            first_shot(self)
            delay(1 * ms)

        if self.t_recooling_after_first_shot_Node1 > 0 or self.t_recooling_after_first_shot_Node2 > 0:
            recooling_after_first_shot(self)

        ### first_shot doesn't turn off the fiber AOMs
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
        if not self.PGC_and_RO_with_on_chip_beams_Node1:
            self.dds_AOM_A5_Node1.sw.off()
            self.dds_AOM_A6_Node1.sw.off()
        if not self.PGC_and_RO_with_on_chip_beams_Node2:
            self.dds_AOM_A5_Node2.sw.off()
            self.dds_AOM_A6_Node2.sw.off()
        delay(5 * us)

        ### The excitation tone on both GRIN2 AOMs, behind closed gates. An
        ### override of p_GRIN2_excitation_NodeX to -150 leaves that node dark.
        ### GRIN1 is set and switched by the pumping trace itself.
        self.ttl_GRIN2_switch_Node1.on()  # turns off excitation
        self.ttl_GRIN2_switch_Node2.on()
        delay(5 * us)
        ### float(): an override or scan value written as -150 arrives as an int
        self.dds_GRIN2_Node1.set(frequency=self.f_GRIN2_excitation_Node1,
                                 amplitude=dB_to_V(float(self.p_GRIN2_excitation_Node1)))
        self.dds_GRIN2_Node2.set(frequency=self.f_GRIN2_excitation_Node2,
                                 amplitude=dB_to_V(float(self.p_GRIN2_excitation_Node2)))
        delay(5 * us)
        self.dds_GRIN2_Node1.sw.on()
        self.dds_GRIN2_Node2.sw.on()
        delay(5 * us)

        excitation_cycle = 0
        for excitation_cycle in range(max_excitation_cycles):
            self.core.break_realtime()

            ############################### CW optical pumping, both nodes: F=1,m_F=0
            if do_pumping:
                self.core_dma.playback_handle(OP_handle)

            ############################### excitation F=1,m=0 -> F'=0,m'=0; photons on the shared SPCMs
            for excitation_attempt in range(n_excitation_attempts):
                ### timestamp_mu below waits for the gate to close, so each
                ### attempt starts with no slack
                slack = now_mu() - self.core.get_rtio_counter_mu()
                if slack < 1e5:
                    self.core.break_realtime()

                dma_start_mu = now_mu()
                self.core_dma.playback_handle(excitation_handle)

                t1 = dma_start_mu + self.core.seconds_to_mu(5 * us)
                t_end_SPCM = t1 + gate_end_mu

                SPCM0_click_time = self.ttl_SPCM0.timestamp_mu(t_end_SPCM)
                SPCM1_click_time = self.ttl_SPCM1.timestamp_mu(t_end_SPCM)
                SPCM0_OtherNode_click_time = self.ttl_SPCM0_OtherNode.timestamp_mu(t_end_SPCM)
                SPCM1_OtherNode_click_time = self.ttl_SPCM1_OtherNode.timestamp_mu(t_end_SPCM)

                if (SPCM0_click_time > 0 or SPCM1_click_time > 0
                        or SPCM0_OtherNode_click_time > 0 or SPCM1_OtherNode_click_time > 0):
                    if Any_SPCM_click_counter < short_timestamps_length:
                        if SPCM0_click_time > 0:
                            SPCM0_timestamps[Any_SPCM_click_counter] = self.core.mu_to_seconds(SPCM0_click_time)
                        if SPCM1_click_time > 0:
                            SPCM1_timestamps[Any_SPCM_click_counter] = self.core.mu_to_seconds(SPCM1_click_time)
                        if SPCM0_OtherNode_click_time > 0:
                            SPCM0_OtherNode_timestamps[Any_SPCM_click_counter] = self.core.mu_to_seconds(SPCM0_OtherNode_click_time)
                        if SPCM1_OtherNode_click_time > 0:
                            SPCM1_OtherNode_timestamps[Any_SPCM_click_counter] = self.core.mu_to_seconds(SPCM1_OtherNode_click_time)
                        tStamps_t1[Any_SPCM_click_counter] = self.core.mu_to_seconds(t1)
                        Any_SPCM_click_counter += 1

                    ### Drain all additional edges, including when the row buffer is already full
                    if SPCM0_click_time > 0: self.ttl_SPCM0.count(t_end_SPCM)
                    if SPCM1_click_time > 0: self.ttl_SPCM1.count(t_end_SPCM)
                    if SPCM0_OtherNode_click_time > 0: self.ttl_SPCM0_OtherNode.count(t_end_SPCM)
                    if SPCM1_OtherNode_click_time > 0: self.ttl_SPCM1_OtherNode.count(t_end_SPCM)

                delay(5 * us)

            ############################ recooling, both nodes
            if do_recooling and (excitation_cycle + 1) % recool_every_n_OP == 0:
                self.core_dma.playback_handle(recooling_handle)
                delay(5 * us)

                # ### Back to the OP coils early, as single_photon_experiment_9
                # ### does: the pumping trace then finds the field already
                # ### moving from PGC to OP, the condition single-node pumping
                # ### was characterised under.
                # self.zotino0_Node1.set_dac(
                #     [self.AZ_bottom_volts_OP_Node1, -self.AZ_bottom_volts_OP_Node1, self.AX_volts_OP_Node1, self.AY_volts_OP_Node1],
                #     channels=self.coil_channels_Node1)
                # self.zotino0_Node2.set_dac(
                #     [self.AZ_bottom_volts_OP_Node2, -self.AZ_bottom_volts_OP_Node2, self.AX_volts_OP_Node2, self.AY_volts_OP_Node2],
                #     channels=self.coil_channels_Node2)
                # delay(0.4 * ms)

            ############################# both atoms still there? every atom_check_every_n cycles
            if (excitation_cycle + 1) % atom_check_every_n == 0:
                ### set_dac writes ~4.7 us in the past, and with no recooling
                ### before it this follows the last photon read with almost no
                ### slack. The check has no timing relation to the excitation.
                self.core.break_realtime()
                self.zotino0_Node1.set_dac(
                    [self.AZ_bottom_volts_PGC_Node1, -self.AZ_bottom_volts_PGC_Node1, self.AX_volts_PGC_Node1, self.AY_volts_PGC_Node1],
                    channels=self.coil_channels_Node1)
                self.zotino0_Node2.set_dac(
                    [self.AZ_bottom_volts_PGC_Node2, -self.AZ_bottom_volts_PGC_Node2, self.AX_volts_PGC_Node2, self.AY_volts_PGC_Node2],
                    channels=self.coil_channels_Node2)
                delay(0.4 * ms)

                self.dds_cooling_DP_Node1.set(frequency=self.f_cooling_DP_RO_Node1, amplitude=self.ampl_cooling_DP_RO_Node1)
                self.dds_cooling_DP_Node2.set(frequency=self.f_cooling_DP_RO_Node2, amplitude=self.ampl_cooling_DP_RO_Node2)
                delay(5 * us)

                ### fiber AOMs first, then the light reaches both atoms on one timestamp
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
                delay(5 * us)
                self.dds_cooling_DP_Node1.sw.on()
                self.dds_cooling_DP_Node2.sw.on()
                self.ttl_repump_switch_Node1.off()
                self.ttl_repump_switch_Node2.off()
                delay(1 * us)

                with parallel:
                    self.ttl_SPCM0_counter.gate_rising(t_atom_check)
                    self.ttl_SPCM1_counter.gate_rising(t_atom_check)
                    self.ttl_SPCM0_OtherNode_counter.gate_rising(t_atom_check)
                    self.ttl_SPCM1_OtherNode_counter.gate_rising(t_atom_check)

                AllSPCMs_RO_atom_check = int(self.ttl_SPCM0_counter.fetch_count()
                                             + self.ttl_SPCM1_counter.fetch_count()
                                             + self.ttl_SPCM0_OtherNode_counter.fetch_count()
                                             + self.ttl_SPCM1_OtherNode_counter.fetch_count())
                AllSPCMs_RO_atom_check_array[int(excitation_cycle / atom_check_every_n)] = AllSPCMs_RO_atom_check

                ### stop the excitation once either atom is lost
                if AllSPCMs_RO_atom_check / t_atom_check < self.two_atom_threshold:
                    delay(100 * us)
                    break

                delay(10 * us)
                self.dds_cooling_DP_Node1.sw.off()
                self.dds_cooling_DP_Node2.sw.off()
                self.ttl_repump_switch_Node1.on()
                self.ttl_repump_switch_Node2.on()
                delay(5 * us)
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
                if not self.PGC_and_RO_with_on_chip_beams_Node1:
                    self.dds_AOM_A5_Node1.sw.off()
                    self.dds_AOM_A6_Node1.sw.off()
                if not self.PGC_and_RO_with_on_chip_beams_Node2:
                    self.dds_AOM_A5_Node2.sw.off()
                    self.dds_AOM_A6_Node2.sw.off()

                ### back to the OP coils, as single_photon_experiment_9 does
                ### after a passed check (the gate read above left no slack)
                self.core.break_realtime()
                # self.zotino0_Node1.set_dac(
                #     [self.AZ_bottom_volts_OP_Node1, -self.AZ_bottom_volts_OP_Node1, self.AX_volts_OP_Node1, self.AY_volts_OP_Node1],
                #     channels=self.coil_channels_Node1)
                # self.zotino0_Node2.set_dac(
                #     [self.AZ_bottom_volts_OP_Node2, -self.AZ_bottom_volts_OP_Node2, self.AX_volts_OP_Node2, self.AY_volts_OP_Node2],
                #     channels=self.coil_channels_Node2)
                # delay(0.4 * ms)

            delay(5 * us)

        ### The traces baked in the A5/A6 amplitudes at record time and the
        ### loader may have run feedback since: restore the live ones while
        ### the RF switches are off. Resynced first: four DDS writes per node
        ### would otherwise run on the few us left after the last photon read.
        self.core.break_realtime()
        self.dds_AOM_A5_Node1.set(frequency=self.AOM_A5_freq_Node1, amplitude=self.stabilizer_AOM_A5_Node1.amplitude)
        self.dds_AOM_A6_Node1.set(frequency=self.AOM_A6_freq_Node1, amplitude=self.stabilizer_AOM_A6_Node1.amplitude)
        self.dds_AOM_A5_Node2.set(frequency=self.AOM_A5_freq_Node2, amplitude=self.stabilizer_AOM_A5_Node2.amplitude)
        self.dds_AOM_A6_Node2.set(frequency=self.AOM_A6_freq_Node2, amplitude=self.stabilizer_AOM_A6_Node2.amplitude)
        delay(5 * us)

        self.dds_GRIN2_Node1.sw.off()
        self.dds_GRIN2_Node2.sw.off()
        delay(0.1 * ms)

        second_shot(self)

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
        self.dds_AOM_A5_Node1.sw.off()
        self.dds_AOM_A5_Node2.sw.off()
        self.dds_AOM_A6_Node1.sw.off()
        self.dds_AOM_A6_Node2.sw.off()
        delay(0.1 * ms)

        end_measurement(self)

        self.core.break_realtime()

        ### only the passed checks: the loop exits on the failed one
        for val in AllSPCMs_RO_atom_check_array[0:int(excitation_cycle / atom_check_every_n)]:
            self.append_to_dataset('AllSPCMs_RO_atom_check', val)

        for i in range(Any_SPCM_click_counter):
            self.append_to_dataset('SPCM0_SinglePhoton_reduced_tStamps', SPCM0_timestamps[i])
            self.append_to_dataset('SPCM1_SinglePhoton_reduced_tStamps', SPCM1_timestamps[i])
            self.append_to_dataset('SPCM0_OtherNode_SinglePhoton_reduced_tStamps', SPCM0_OtherNode_timestamps[i])
            self.append_to_dataset('SPCM1_OtherNode_SinglePhoton_reduced_tStamps', SPCM1_OtherNode_timestamps[i])
            self.append_to_dataset('reference_tStamps_t1', tStamps_t1[i])

        ### Number of compact timestamp rows saved during this measurement
        self.append_to_dataset('n_photon_events', Any_SPCM_click_counter)

        self.append_to_dataset('n_excitation_cycles', excitation_cycle)
        self.core.break_realtime()

    self.append_to_dataset('n_feedback_per_iteration', self.n_feedback_per_iteration)
    self.append_to_dataset('n_atom_loaded_per_iteration', self.n_atom_loaded_per_iteration)
