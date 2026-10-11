"""Measure the deterministic master-satellite TTL skew, in both directions.

WHY THIS EXISTS ALONGSIDE measure_drtio_event_cost.py
-----------------------------------------------------
That probe (RID 38429, 2026-09-10) answered "will an event be accepted":
submission cost 713 vs 711 mu/event and an identical 1000 mu minimum viable
lead on both destinations, because ARTIQ buffers events at the satellite ahead
of time. It is a cabling-free control and must stay that way, so this is a
separate file rather than an extension of it.

It did NOT answer "when does the pin actually move", and it touched only an
OUTPUT channel (set_config_mu -> one rtio_output). Remote INPUT latency has
never been measured in this repository, while plan_codex_detail.md:569 and
PLAN.md:945-952 both flag remote input/TTL-timestamp operations as
latency-sensitive. This probe measures both directions.

WHAT IS BEING MEASURED, AND WHY IT IS TWO DIFFERENT THINGS
----------------------------------------------------------
Cross-node timing has two independent layers, and the pre-DRTIO code conflated
them. This probe measures LAYER 1 only:

  LAYER 1, electrical.  Do the two crates' output PINS move together for one
      scheduled timestamp?  -> t_Node2_rtio_offset_mu.  Everything with a gate
      wider than a microsecond (atom loading, the alternating readout) cares
      only about this.  THIS FILE.

  LAYER 2, optical.  Do the two emitted PHOTONS arrive at the SPCM together?
      -> t_Node2_excitation_delay_mu.  Firing both excitation pulses at the
      same now_mu does not make the light reach the two atoms together: fibre
      and free-space path lengths differ, AOM turn-on differs, and each atom's
      emission is referenced to its own local pulse.  NOT measurable here --
      it needs a photon coincidence measurement (HOM dip / cross-correlation).

CABLING
-------
Both cross-node coax ALREADY EXIST -- they are the legacy TTL handshake
cables, named in BaseExperiment_master_satellite.NODE_TTL_ALIASES. Add two BNC
T-adapters and two short coax:

    role          signal                     from              to
    existing      master out -> satellite    ttl12  (dest 0)   ttl24 (dest 1)
    NEW (T)       master out -> master       ttl12  (dest 0)   ttl3  (dest 0)
    existing      satellite out -> master    ttl31  (dest 1)   ttl11 (dest 0)
    NEW (T)       satellite out -> satellite ttl31  (dest 1)   ttl27 (dest 1)

ttl2/ttl3/ttl10 are free counter-capable TTLInOut on the master; ttl27 is free
on the satellite (Node2's ttl11, "not being used").

The T-adapter is the whole trick: ONE output edge produces TWO timestamps, so
the output path cancels exactly and the difference is pure input-side skew.

Cable length differences between the three cross-node lines are < 3 ns
(measured by the lab), i.e. under half a 125 MHz RTIO clock period, so there
is no need to swap cables and re-run. Recorded here so nobody re-derives it.

WHAT IT REPORTS
---------------
  T1  input skew      pulse ttl12; d_sat - d_mas over the two timestamps.
                      The number that has never existed in this repository.
  T2  output skew     pulse ttl31; same construction the other way. Combined
                      with T1 this isolates how much later a satellite output
                      pin moves for one scheduled timestamp -- the value for
                      t_Node2_rtio_offset_mu.
  T3  consistency     the two round trips must agree within jitter. A
                      mismatch means a cable or channel assumption is wrong,
                      not that DRTIO is bad.
  T4  counter aperture  the loader uses EDGE COUNTERS, not timestamps. Sweep
                      a fixed gate across the pulse in 4 ns steps and find the
                      0->1 count transition on a master and a satellite
                      counter. It must agree with T1, which is what proves the
                      timestamp-measured quantity governs a gated count.

STABILITY IS WHAT DECIDES THE DESIGN
------------------------------------
Deterministic within a link session is not the same as constant. Levels 1-3
are automated here; level 4 is an operator procedure:

  1. repeats inside one kernel invocation          (n_shots)
  2. repeats across kernel invocations in one run  (n_blocks)
  3. across core.reset()                           (n_resets)
  4. ACROSS A SATELLITE LINK RE-ESTABLISH -- reboot the satellite or unplug
     and replug the fibre, then run this again and compare. NOT automated,
     because it needs physical intervention.

Level 4 is the risk. If 1-3 give sub-nanosecond spread but 4 moves by a clock
period (8 ns), the offset is a per-session quantity, not a constant to bake
into a delay(). For atom loading (10 us - 1 ms gates) an 8 ns session shift is
irrelevant and the answer is "insensitive, proceed". For the deferred
herald/photon work it is decisive, which is why both are measured at once.

ACCEPTANCE GATE
---------------
|T1| and |T2| each < 1 us, and the level-1..3 std < 100 ns. Both are ~100x
looser than the physics needs and ~1000x looser than the buffer-ahead
behaviour predicts, so a failure means something structural (clocking, link,
SED lanes), and the two-node loader must not be written until it is
understood.

RESULT (fill in after the first hardware run, as
measure_drtio_event_cost.py:50-73 does -- that habit is why the DRTIO
question is settled today)

    RESULT, <date> (RID <n>)
        T1 input skew      <median> mu   std <...>   n <...>
        T2 output skew     <median> mu   std <...>   n <...>
        T4 aperture        master <...>  satellite <...>
        level 4            <moved / did not move> across a link re-establish
        -> t_Node2_rtio_offset_mu = <value>

    Copy T2's median into t_Node2_rtio_offset_mu WITH ITS SIGN UNCHANGED.
    Positive means Node2's pin moves later for the same timestamp, and the
    sequence cancels that by scheduling Node2 earlier -- the convention is
    spelled out at the variable's declaration in
    ExperimentVariables_master_satellite_global.py. Negating it here would
    double the skew rather than remove it, and the only symptom would be
    readout light leaking between the alternating windows.

SAFETY
------
Drives two TTL outputs that exist solely for the retired inter-node handshake
and reads four TTL inputs. No DDS, no Zotino, no coils, no DMA. Does not use
BaseExperiment, so it runs regardless of the two-node base-layer work. The
core reset in the stability sweep drops TTL output states; that is why it is
behind an argument.
"""

import numpy as np

from artiq.experiment import *


#: Destination 1 readiness, same bounded-poll shape BaseExperiment uses.
_SATELLITE_POLL_ATTEMPTS = 100
_SATELLITE_POLL_INTERVAL_MU = 100_000_000  # 100 ms in mu at 1 ns/mu

#: A missed edge: TTLInOut.timestamp_mu returns -1 when the gate saw nothing.
_NO_EDGE = np.int64(-1)


class measure_drtio_ttl_skew(EnvExperiment):
    """measure_drtio_ttl_skew"""

    def build(self):
        self.setattr_device("core")

        # Outputs. Both are TTLOut and both exist only for the retired
        # handshake, so driving them disturbs nothing.
        self.setattr_device("ttl12")   # dest 0, master   -> ttl24 and ttl3
        self.setattr_device("ttl31")   # dest 1, satellite-> ttl11 and ttl27

        # Inputs. All four are TTLInOut and must be put in input mode.
        self.setattr_device("ttl3")    # dest 0, T off ttl12   (NEW cable)
        self.setattr_device("ttl24")   # dest 1, from ttl12    (existing)
        self.setattr_device("ttl11")   # dest 0, from ttl31    (existing)
        self.setattr_device("ttl27")   # dest 1, T off ttl31   (NEW cable)

        # Edge counters for T4. Only the two that see the master-driven pulse.
        self.setattr_device("ttl3_counter")    # dest 0
        self.setattr_device("ttl24_counter")   # dest 1

        self.setattr_argument(
            "n_shots",
            NumberValue(1000, ndecimals=0, step=100, min=1, max=100000),
            "Stability",
            tooltip="Level 1: repeats inside one kernel invocation.",
        )
        self.setattr_argument(
            "n_blocks",
            NumberValue(10, ndecimals=0, step=1, min=1, max=1000),
            "Stability",
            tooltip="Level 2: kernel invocations per run.",
        )
        self.setattr_argument(
            "n_resets",
            NumberValue(0, ndecimals=0, step=1, min=0, max=50),
            "Stability",
            tooltip="Level 3: extra passes separated by core.reset(). "
                    "WARNING: a core reset drops TTL output states. Leave 0 "
                    "unless levels 1-2 look clean and you want level 3.",
        )
        self.setattr_argument(
            "run_aperture_sweep",
            BooleanValue(True),
            "Counter aperture",
            tooltip="T4. Sweeps a fixed gate across the pulse and finds the "
                    "0->1 count transition on a master and a satellite "
                    "counter. This is the one that speaks to the loader, "
                    "which gates counters rather than reading timestamps.",
        )
        self.setattr_argument(
            "aperture_half_range_ns",
            NumberValue(200, ndecimals=0, step=50, min=20, max=5000),
            "Counter aperture",
        )
        self.setattr_argument(
            "aperture_step_ns",
            NumberValue(4, ndecimals=0, step=1, min=1, max=100),
            "Counter aperture",
            tooltip="4 ns is half an RTIO coarse period at 125 MHz.",
        )
        self.setattr_argument(
            "aperture_gate_ns",
            NumberValue(200, ndecimals=0, step=50, min=20, max=5000),
            "Counter aperture",
        )
        self.setattr_argument(
            "aperture_repeats",
            NumberValue(64, ndecimals=0, step=8, min=1, max=1000),
            "Counter aperture",
        )

    def prepare(self):
        self.n_shots = int(self.n_shots)
        self.n_blocks = int(self.n_blocks)
        self.n_resets = int(self.n_resets)
        self.aperture_repeats = int(self.aperture_repeats)

        # Gate wide enough to contain the pulse comfortably even if the skew
        # turns out to be far larger than expected; the measurement is the
        # timestamp, not the gate, so a wide gate costs only timeline.
        self.gate_mu = np.int64(10_000)        # 10 us
        self.pulse_offset_mu = np.int64(5_000)  # fire mid-gate
        self.pulse_mu = np.int64(1_000)         # 1 us, well above any jitter

        # Per-shot results, filled in the kernel. int64 because they are
        # machine units and may be negative.
        self.d_master = np.full(self.n_shots, _NO_EDGE, dtype=np.int64)
        self.d_satellite = np.full(self.n_shots, _NO_EDGE, dtype=np.int64)

        self.aperture_offsets = [
            int(value)
            for value in range(
                -int(self.aperture_half_range_ns),
                int(self.aperture_half_range_ns) + 1,
                int(self.aperture_step_ns),
            )
        ]
        self.aperture_gate_mu = np.int64(int(self.aperture_gate_ns))

        # Collected across blocks/resets, on the host.
        self.t1_samples = []   # input skew, satellite minus master
        self.t2_samples = []   # output skew round trip
        self.miss_count = 0
        self.shot_count = 0

    # ------------------------------------------------------------------ RPC

    def collect(self, which, n_valid):
        """Pull one kernel block's results onto the host.

        Called from the kernel after each block so the arrays can be reused.
        `which` is 1 for the master-driven direction and 2 for the
        satellite-driven one.
        """
        samples = []
        for index in range(int(n_valid)):
            master = int(self.d_master[index])
            satellite = int(self.d_satellite[index])
            if master == -1 or satellite == -1:
                self.miss_count += 1
                continue
            samples.append(satellite - master)
        self.shot_count += int(n_valid)
        if which == 1:
            self.t1_samples += samples
        else:
            self.t2_samples += samples

    def record_aperture(self, tag, offset_ns, count):
        self.aperture_results.setdefault(int(tag), []).append(
            (int(offset_ns), int(count))
        )

    # --------------------------------------------------------------- kernels

    @kernel
    def prepare_hardware(self):
        self.core.reset()

        # Destination 1 must be up before any satellite channel is touched.
        ready = False
        for _ in range(_SATELLITE_POLL_ATTEMPTS):
            if self.core.get_rtio_destination_status(1):
                ready = True
                break
            delay_mu(_SATELLITE_POLL_INTERVAL_MU)
        if not ready:
            raise RuntimeError("DRTIO destination 1 did not become available.")
        self.core.break_realtime()

        # Inputs in input mode; outputs driven low so the first edge of a shot
        # is unambiguous.
        self.ttl3.input()
        self.ttl24.input()
        self.ttl11.input()
        self.ttl27.input()
        delay(1 * ms)
        self.ttl12.off()
        self.ttl31.off()
        delay(1 * ms)

    @kernel
    def shots_master_drives(self, n: TInt32):
        """Pulse the MASTER output; timestamp on a master and a satellite in.

        d_master     = ts(ttl3)  - t_fire   master out -> master in
        d_satellite  = ts(ttl24) - t_fire   master out -> satellite in
        T1 = d_satellite - d_master, so the output path cancels and what is
        left is the input-side skew plus the (<3 ns) cable difference.
        """
        for shot in range(n):
            self.core.break_realtime()
            t_fire = np.int64(0)
            with parallel:
                gate_master = self.ttl3.gate_rising_mu(self.gate_mu)
                gate_satellite = self.ttl24.gate_rising_mu(self.gate_mu)
                with sequential:
                    delay_mu(self.pulse_offset_mu)
                    t_fire = now_mu()
                    self.ttl12.pulse_mu(self.pulse_mu)

            stamp_master = self.ttl3.timestamp_mu(gate_master)
            stamp_satellite = self.ttl24.timestamp_mu(gate_satellite)

            if stamp_master == _NO_EDGE:
                self.d_master[shot] = _NO_EDGE
            else:
                self.d_master[shot] = stamp_master - t_fire
            if stamp_satellite == _NO_EDGE:
                self.d_satellite[shot] = _NO_EDGE
            else:
                self.d_satellite[shot] = stamp_satellite - t_fire

    @kernel
    def shots_satellite_drives(self, n: TInt32):
        """Pulse the SATELLITE output; timestamp on a satellite and a master in.

        d_master    = ts(ttl11) - t_fire   satellite out -> master in
        d_satellite = ts(ttl27) - t_fire   satellite out -> satellite in
        Compared with T1 this says how much later the satellite's output pin
        moves for the same scheduled timestamp.
        """
        for shot in range(n):
            self.core.break_realtime()
            t_fire = np.int64(0)
            with parallel:
                gate_master = self.ttl11.gate_rising_mu(self.gate_mu)
                gate_satellite = self.ttl27.gate_rising_mu(self.gate_mu)
                with sequential:
                    delay_mu(self.pulse_offset_mu)
                    t_fire = now_mu()
                    self.ttl31.pulse_mu(self.pulse_mu)

            stamp_master = self.ttl11.timestamp_mu(gate_master)
            stamp_satellite = self.ttl27.timestamp_mu(gate_satellite)

            if stamp_master == _NO_EDGE:
                self.d_master[shot] = _NO_EDGE
            else:
                self.d_master[shot] = stamp_master - t_fire
            if stamp_satellite == _NO_EDGE:
                self.d_satellite[shot] = _NO_EDGE
            else:
                self.d_satellite[shot] = stamp_satellite - t_fire

    @kernel
    def aperture_point(self, offset_ns: TInt32, repeats: TInt32):
        """One rung of T4: fixed gate placed offset_ns from the pulse.

        Both counters are gated identically relative to the SAME master-driven
        pulse, so the difference between the two transition positions is the
        same quantity T1 measures -- but measured the way the loader actually
        works, through EdgeCounter rather than timestamp_mu.
        """
        self.core.break_realtime()
        offset_mu = np.int64(offset_ns)

        for _ in range(repeats):
            self.core.break_realtime()
            with parallel:
                with sequential:
                    delay_mu(self.pulse_offset_mu + offset_mu)
                    self.ttl3_counter.gate_rising_mu(self.aperture_gate_mu)
                with sequential:
                    delay_mu(self.pulse_offset_mu + offset_mu)
                    self.ttl24_counter.gate_rising_mu(self.aperture_gate_mu)
                with sequential:
                    delay_mu(self.pulse_offset_mu)
                    self.ttl12.pulse_mu(self.pulse_mu)

        count_master = 0
        count_satellite = 0
        for _ in range(repeats):
            count_master += self.ttl3_counter.fetch_count()
            count_satellite += self.ttl24_counter.fetch_count()

        self.record_aperture(0, offset_ns, count_master)
        self.record_aperture(1, offset_ns, count_satellite)

    # ------------------------------------------------------------------ host

    def run(self):
        self.aperture_results = {}

        print("measure_drtio_ttl_skew")
        print("  cable deltas between the cross-node lines are < 3 ns "
              "(lab measurement), i.e. under half an RTIO coarse period")
        print("  %d shots x %d blocks, %d extra core-reset passes"
              % (self.n_shots, self.n_blocks, self.n_resets))

        self.prepare_hardware()

        for pass_index in range(self.n_resets + 1):
            if pass_index > 0:
                print("  -- core reset, pass %d" % pass_index)
                self.prepare_hardware()
            for _ in range(self.n_blocks):
                # The kernels fill d_master/d_satellite in place and
                # return nothing; the host already knows the shot count, so
                # there is no reason to round-trip it through a return value.
                self.shots_master_drives(self.n_shots)
                self.collect(1, self.n_shots)
                self.shots_satellite_drives(self.n_shots)
                self.collect(2, self.n_shots)

        if self.run_aperture_sweep:
            print("  T4 aperture sweep: %d points x %d repeats"
                  % (len(self.aperture_offsets), self.aperture_repeats))
            for offset_ns in self.aperture_offsets:
                self.aperture_point(offset_ns, self.aperture_repeats)

        self.report()

    def _stats(self, samples):
        if not samples:
            return None
        array = np.array(samples, dtype=np.float64)
        return {
            "n": int(array.size),
            "median": float(np.median(array)),
            "mean": float(array.mean()),
            "std": float(array.std()),
            "min": float(array.min()),
            "max": float(array.max()),
        }

    def _transition_ns(self, tag):
        """Where the gate first stops catching the pulse, in ns."""
        points = sorted(self.aperture_results.get(tag, []))
        if not points:
            return None
        peak = max(count for _, count in points)
        if peak == 0:
            return None
        half = peak / 2.0
        for (offset, count), (next_offset, next_count) in zip(points, points[1:]):
            if count >= half > next_count:
                return float(offset + next_offset) / 2.0
        return None

    def report(self):
        print("")
        print("=" * 68)
        t1 = self._stats(self.t1_samples)
        t2 = self._stats(self.t2_samples)

        for name, stats in (("T1 input skew  (sat-in minus mas-in)", t1),
                            ("T2 output skew (sat-out round trip) ", t2)):
            if stats is None:
                print("  %s   NO VALID SHOTS" % name)
                continue
            print("  %s" % name)
            print("      n=%d  median=%.1f mu  mean=%.1f  std=%.1f  "
                  "min=%.0f  max=%.0f"
                  % (stats["n"], stats["median"], stats["mean"],
                     stats["std"], stats["min"], stats["max"]))

        if self.miss_count:
            print("  MISSED EDGES: %d of %d shots had a -1 timestamp. A gate "
                  "missed the pulse; widen gate_mu or check the cabling."
                  % (self.miss_count, self.shot_count))

        if t1 is not None and t2 is not None:
            # T2 - T1 removes the common input-side term, leaving the
            # output-side asymmetry: the candidate for t_Node2_rtio_offset_mu.
            offset = t2["median"] - t1["median"]
            print("")
            print("  => t_Node2_rtio_offset_mu candidate: %.1f mu" % offset)
            print("     (T2 median minus T1 median; the shared input-side "
                   "term cancels)")

        if self.run_aperture_sweep:
            master_edge = self._transition_ns(0)
            satellite_edge = self._transition_ns(1)
            print("")
            print("  T4 counter aperture (the quantity the loader gates on)")
            if master_edge is None or satellite_edge is None:
                print("      transition not found; check the sweep range")
            else:
                print("      master    transition at %+.1f ns" % master_edge)
                print("      satellite transition at %+.1f ns" % satellite_edge)
                print("      difference %+.1f ns  -- must agree with T1"
                      % (satellite_edge - master_edge))

        # Acceptance gate, stated rather than enforced: a failure here means
        # something structural, and the operator should see the numbers.
        print("")
        if t1 is not None and t2 is not None:
            within = (abs(t1["median"]) < 1000 and abs(t2["median"]) < 1000
                      and t1["std"] < 100 and t2["std"] < 100)
            print("  ACCEPTANCE (|skew| < 1000 mu, std < 100 mu): %s"
                  % ("PASS" if within else "FAIL -- investigate before "
                     "writing two-node sequences"))
        print("  LEVEL 4 is not automated: reboot the satellite or replug the "
              "fibre, re-run, and compare the medians above.")
        print("=" * 68)

        for name, stats in (("input", t1), ("output", t2)):
            if stats is None:
                continue
            self.set_dataset("drtio_skew_%s_mu" % name,
                             stats["median"], broadcast=True, persist=True)
            self.set_dataset("drtio_skew_%s_std_mu" % name,
                             stats["std"], broadcast=True, persist=True)
            self.set_dataset("drtio_skew_%s_n" % name,
                             stats["n"], broadcast=True, persist=True)
        if self.run_aperture_sweep:
            master_edge = self._transition_ns(0)
            satellite_edge = self._transition_ns(1)
            if master_edge is not None and satellite_edge is not None:
                self.set_dataset("drtio_skew_counter_aperture_mu",
                                 satellite_edge - master_edge,
                                 broadcast=True, persist=True)
