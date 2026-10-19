"""measure_drtio_ttl_skew.py

Print the measured TTL delay of every cabled path between and within the two
nodes, averaged over n_shots (default 1000). Nothing is estimated, derived or
corrected: every printed number is a plain statistic of measured timestamps.

WHAT ONE NUMBER MEANS
---------------------
    delay = (input timestamp) - (scheduled fire timestamp of the output)

in machine units; 1 mu = 1 ns (ref_period 1 ns). The fire timestamp is the
timeline position of the output's rising edge (now_mu() at the pulse). The
input timestamp is the first rising edge that input saw inside that output's
window. The DIO TTL channels resolve 1 ns on both inputs and outputs.

Each timestamp is in the clock of the node that recorded it, and a Node2
output fires when NODE2's clock reaches its timestamp. Legs 4 and 5 stay
inside one node. Legs 1-3 compare the two nodes' clocks, so they include the
offset between those clocks. Context, read from the ARTIQ 7 DRTIO gateware:
Node1 sends its raw timestamp counter when the link comes up
(gateware/drtio/rt_controller_master.py:54-55) and Node2 loads it as received
(rt_packet_satellite.py:85-86), with no compensation for the link latency,
so Node2's clock reads behind Node1's. Only a link re-establish reloads it
(artiq-zynq release-7 src/runtime/src/rtio_mgt.rs:276-286); core.reset()
does not. So a negative delay on legs 1-2 is not an error, and cross-node
numbers belong to the link session they were measured in. To see whether
legs 1-3 repeat from one link session to the next, reboot Node2 or replug
the fibre, run this again and compare the printed means.

LEGS AND CABLING
----------------
Unified device names; Node2's local ttlN is unified ttl(N+16).

    leg  direction        fires  stamps  physical connection
    1    Node1 -> Node2   ttl12  ttl25   cable 1: Node1 ttl12 -> Node2 ttl9
    2    Node1 -> Node2   ttl15  ttl26   cable 2: Node1 ttl15 -> Node2 ttl10
    3    Node2 -> Node1   ttl31  ttl11   cable 3: Node2 ttl15 -> Node1 ttl11
    4    Node1 -> Node1   ttl12  ttl10   loopback (add it, see below)
    5    Node2 -> Node2   ttl31  ttl27   loopback (add it, see below)

Cables 1-3 exist, wired as the lab stated on 2026-10-18. Two repo records
put cables 1 and 2 one Node2 input lower, on Node2 ttl8 and ttl9 (ttl24,
ttl25): NODE_TTL_ALIASES (utilities/BaseExperiment_master_satellite.py:
181-182) and the test_ttl_pulse_experiment docstring
(subroutines/experiment_functions.py:88-90). If leg 1 or 2 sees no edge,
check which Node2 input that cable is in.

Each loopback is a BNC T-adapter on the output plus a short coax into a free
input on the same node:

    leg 4   T on Node1 ttl12, short coax into Node1 ttl10
    leg 5   T on Node2 ttl15 (ttl31), short coax into Node2 ttl11 (ttl27)

Legs 1 and 4 see the same ttl12 edge; legs 3 and 5 see the same ttl31 edge.
Without the loopbacks, legs 4 and 5 print "no edge" and legs 1-3 still run;
an edge on leg 4 or 5 without its loopback means something else is plugged
into that input (no current code uses Node1 ttl10 or Node2 ttl11).
The Node1 loop uses ttl10, not ttl3: ttl3 is the GRIN1 toggle-switch input
read by SamplerMOTCoilAndBeamBalanceTune_master_satellite.py:392-395.

ONE SHOT
--------
Each output fires once, in its own 10 us window, 5 us into the window. Only
the input(s) wired to that output are gated, and only during that window, so
a cable on the wrong input shows up as a missed edge, not as a wrong number.
Delays from -5 us to +5 us are caught. After the first edge is read, the
rest of that gate is drained and counted as extra edges (ringing, e.g. on a
T stub, shows up there). A gate that saw no edge is a missed shot and is
left out of the statistics.

PRINTOUT AND DATASETS
---------------------
Three header lines (name and n_shots, what "delay" means), then per leg:
the mean delay over the shots that saw an edge (ns, 2 decimals),
n_valid/n_shots, std (population, ddof=0), min, max, and the total number of
extra edges. A leg with no edge at all prints one line saying so and naming
the connection it expects.

Archived to this run's h5 only (broadcast=False, persist=False), k = 1..5:

    drtio_ttl_skew_leg<k>_delay_mu      int64[n_shots]  delay of each shot
    drtio_ttl_skew_leg<k>_edge_seen     bool[n_shots]   False = missed shot
    drtio_ttl_skew_leg<k>_extra_edges   int32[n_shots]  edges after the first
    drtio_ttl_skew_leg<k>_mean_mu       float           NaN if no edge at all

A missed shot holds 0 in delay_mu, so always select with edge_seen. -1 is not
used as a marker because a real cross-node delay can be -1.

FAILURES
--------
The arrays reach the host only when run_shots returns (artiq-zynq
src/runtime/src/kernel/core1.rs:121-125), so an exception during the shots
(an input FIFO overflowing past 64 queued edges, the link dropping) loses
that whole run and prints no legs; submit again. An input overflow left
over from an earlier experiment is cleared by the start-up drain.

SAFETY
------
Drives only Node1 ttl12, Node1 ttl15 and Node2 ttl15 (ttl31): set low at the
start, one 1 us pulse each per shot, left low. They belong to the retired
inter-node handshake and are wired only to the inputs above. Puts Node1
ttl10, ttl11 and Node2 ttl9, ttl10, ttl11 (ttl25, ttl26, ttl27) in input mode
and gates them. Nothing else: no DDS, Zotino, Sampler or DMA, and Node1's
SPCM inputs ttl0, ttl1, ttl8, ttl9 (all four SPCMs are on Node1) are not
touched. Standalone EnvExperiment, no BaseExperiment. Writes no broadcast or
persistent dataset.

No core.reset(). On Kasli-SoC it resets only Node1's RTIO core (artiq-zynq
src/runtime/src/rtio_csr.rs:22-26). The routine that would also reset the
satellite (rtio_mgt.rs:356-362, which sends it a ResetRequest at :324) has
no caller, so a reset could not clear Node2's inputs. The probe drains its
five inputs itself at the start.

RESULT, 2026-10-18 (RID 39102, loopbacks fitted)
------------------------------------------------
    measure_drtio_ttl_skew  n_shots 1000
    leg 1  Node1 -> Node2  ttl12 -> ttl25   mean -117.00 ns  1000/1000, std 0
    leg 2  Node1 -> Node2  ttl15 -> ttl26   mean -118.00 ns  1000/1000, std 0
    leg 3  Node2 -> Node1  ttl31 -> ttl11   mean  444.17 ns  1000/1000, std 0.37
    leg 4  Node1 -> Node1  ttl12 -> ttl10   mean  159.78 ns  1000/1000, std 0.41
    leg 5  Node2 -> Node2  ttl31 -> ttl27   mean  161.00 ns  1000/1000, std 0

Leg 3 read 444 or 445 and leg 4 read 159 or 160; every other leg gave one
integer on all 1000 shots. Leg 4 saw exactly one extra rising edge in every
shot; no other leg saw any. The connectivity check just before it (RID 39101)
showed the same thing: 198 rising edges on ttl10 from 100 ttl12 pulses, every
other path 100/100, no crosstalk, nothing in the no-pulse window.

The extra edge is at the END of the ttl12 pulse (edge trace, RID 39103):
ttl10 falls at +1166 ns, rises again at +1175 ns and falls at +1182 ns,
relative to the ttl12 rising edge, in every shot. ttl25 on the other T
branch shows no such bounce. The leg 4 rising edge itself is unaffected.

Earlier the same day, before the BNC T pieces were fitted (RID 39100, cabling
checked by RID 39099): leg 1 -122 ns, leg 2 -119 ns, leg 3 +439 ns, each the
identical integer on all 1000 shots, no extra edges.
"""

import numpy as np

from artiq.experiment import *


#: Node2 is DRTIO destination 1. Wait for it for up to 100 x 100 ms.
_SATELLITE_DESTINATION = 1
_SATELLITE_POLL_ATTEMPTS = 100
_SATELLITE_POLL_INTERVAL_MU = np.int64(100_000_000)

#: One window per output per shot, with the edge fired _FIRE_AT_MU into it,
#: so delays from -_FIRE_AT_MU to _WINDOW_MU - _FIRE_AT_MU are caught.
_WINDOW_MU = np.int64(10_000)
_FIRE_AT_MU = np.int64(5_000)
_PULSE_MU = np.int64(1_000)

#: TTLInOut.timestamp_mu returns -1 when the gate saw no edge.
_NO_EDGE = np.int64(-1)

#: Outputs in window order; row index into fire_mu.
_TTL12 = 0
_TTL15 = 1
_TTL31 = 2
_N_OUTPUTS = 3

#: (leg, direction, output, input, output row, connection). Leg k is row
#: k - 1 of stamp_mu and extra_edges.
_LEGS = (
    (1, "Node1 -> Node2", "ttl12", "ttl25", _TTL12,
     "cable 1: Node1 ttl12 -> Node2 ttl9 (ttl25)"),
    (2, "Node1 -> Node2", "ttl15", "ttl26", _TTL15,
     "cable 2: Node1 ttl15 -> Node2 ttl10 (ttl26)"),
    (3, "Node2 -> Node1", "ttl31", "ttl11", _TTL31,
     "cable 3: Node2 ttl15 (ttl31) -> Node1 ttl11"),
    (4, "Node1 -> Node1", "ttl12", "ttl10", _TTL12,
     "loopback: BNC T on Node1 ttl12, short coax into Node1 ttl10"),
    (5, "Node2 -> Node2", "ttl31", "ttl27", _TTL31,
     "loopback: BNC T on Node2 ttl15 (ttl31), short coax into "
     "Node2 ttl11 (ttl27)"),
)
_N_LEGS = len(_LEGS)


class measure_drtio_ttl_skew(EnvExperiment):
    """measure_drtio_ttl_skew"""

    def build(self):
        self.setattr_device("core")

        # Outputs (TTLOut).
        self.setattr_device("ttl12")   # Node1 ttl12 -> legs 1 and 4
        self.setattr_device("ttl15")   # Node1 ttl15 -> leg 2
        self.setattr_device("ttl31")   # Node2 ttl15 -> legs 3 and 5

        # Inputs (TTLInOut).
        self.setattr_device("ttl25")   # Node2 ttl9,  leg 1
        self.setattr_device("ttl26")   # Node2 ttl10, leg 2
        self.setattr_device("ttl11")   # Node1 ttl11, leg 3
        self.setattr_device("ttl10")   # Node1 ttl10, leg 4 (loopback)
        self.setattr_device("ttl27")   # Node2 ttl11, leg 5 (loopback)

        self.setattr_argument(
            "n_shots",
            NumberValue(1000, ndecimals=0, step=100, min=1, max=10000),
            tooltip="Shots to average. Every shot fires ttl12, ttl15 and "
                    "ttl31 once each.",
        )

    def prepare(self):
        self.n_shots = int(self.n_shots)
        # Filled by run_shots. Flat, one array per kind:
        #   stamp_mu[(leg - 1) * n_shots + shot]     input timestamp, -1 none
        #   extra_edges[(leg - 1) * n_shots + shot]  edges after the first
        #   fire_mu[output * n_shots + shot]         scheduled rising edge
        # Raw timestamps are always positive, so -1 is unambiguous HERE; the
        # delays computed from them on the host use a mask instead.
        self.stamp_mu = np.full(_N_LEGS * self.n_shots, _NO_EDGE,
                                dtype=np.int64)
        self.extra_edges = np.zeros(_N_LEGS * self.n_shots, dtype=np.int32)
        self.fire_mu = np.zeros(_N_OUTPUTS * self.n_shots, dtype=np.int64)

    # --------------------------------------------------------------- kernels

    @kernel
    def prepare_hardware(self):
        self.core.break_realtime()

        # Node2 must be up before any of its channels is touched.
        ready = False
        for _ in range(_SATELLITE_POLL_ATTEMPTS):
            if self.core.get_rtio_destination_status(_SATELLITE_DESTINATION):
                ready = True
                break
            # delay_mu alone only moves the timeline cursor; waiting for the
            # RTIO counter to reach it is what makes each attempt take time.
            delay_mu(_SATELLITE_POLL_INTERVAL_MU)
            self.core.wait_until_mu(now_mu())
        if not ready:
            raise RuntimeError("DRTIO destination 1 (Node2) is not up.")
        self.core.break_realtime()

        self.ttl10.input()
        self.ttl11.input()
        self.ttl25.input()
        self.ttl26.input()
        self.ttl27.input()
        delay(1 * ms)
        self.ttl12.off()
        self.ttl15.off()
        self.ttl31.off()
        delay(1 * ms)

        # No gate is open yet, so this only discards edges an earlier
        # experiment gated on these inputs and never read.
        self._drain(self.ttl10)
        self._drain(self.ttl11)
        self._drain(self.ttl25)
        self._drain(self.ttl26)
        self._drain(self.ttl27)

    @kernel
    def _drain(self, ttl):
        """Discard queued input edges, surviving a stale overflow flag."""
        try:
            ttl.count(now_mu())
        except RTIOOverflow:
            # That read cleared the flag without consuming an event
            # (gateware/rtio/input_collector.py:109-117); read the rest.
            ttl.count(now_mu())

    @kernel
    def _read_gate(self, ttl, gate_end_mu: TInt64, leg: TInt32,
                   shot: TInt32):
        """First rising edge of one gate, then drain and count the rest."""
        slot = (leg - 1) * self.n_shots + shot
        stamp = ttl.timestamp_mu(gate_end_mu)
        self.stamp_mu[slot] = stamp
        if stamp != _NO_EDGE:
            # Left in the FIFO, a second edge would become a later read's
            # "first" edge. A -1 means the gate already closed empty.
            self.extra_edges[slot] = ttl.count(gate_end_mu)

    @kernel
    def run_shots(self):
        n = self.n_shots
        for shot in range(n):
            # Fresh slack every shot: the reads below wait for the gates to
            # close, and a Node2 read waits for NODE2's clock, which lags.
            self.core.break_realtime()

            # At most 3 events share a timestamp (at a window boundary),
            # well under the 8-event limit.

            # Window 1: ttl12. Cable 1 -> ttl25; T loopback -> ttl10.
            with parallel:
                gate_end_1 = self.ttl25.gate_rising_mu(_WINDOW_MU)
                self.ttl10.gate_rising_mu(_WINDOW_MU)
                with sequential:
                    delay_mu(_FIRE_AT_MU)
                    self.fire_mu[_TTL12 * n + shot] = now_mu()
                    self.ttl12.pulse_mu(_PULSE_MU)

            # Window 2: ttl15. Cable 2 -> ttl26.
            with parallel:
                gate_end_2 = self.ttl26.gate_rising_mu(_WINDOW_MU)
                with sequential:
                    delay_mu(_FIRE_AT_MU)
                    self.fire_mu[_TTL15 * n + shot] = now_mu()
                    self.ttl15.pulse_mu(_PULSE_MU)

            # Window 3: ttl31 on Node2. Cable 3 -> ttl11; T loopback -> ttl27.
            with parallel:
                gate_end_3 = self.ttl11.gate_rising_mu(_WINDOW_MU)
                self.ttl27.gate_rising_mu(_WINDOW_MU)
                with sequential:
                    delay_mu(_FIRE_AT_MU)
                    self.fire_mu[_TTL31 * n + shot] = now_mu()
                    self.ttl31.pulse_mu(_PULSE_MU)

            # Read only after every event of the shot is submitted.
            self._read_gate(self.ttl25, gate_end_1, 1, shot)
            self._read_gate(self.ttl10, gate_end_1, 4, shot)
            self._read_gate(self.ttl26, gate_end_2, 2, shot)
            self._read_gate(self.ttl11, gate_end_3, 3, shot)
            self._read_gate(self.ttl27, gate_end_3, 5, shot)

    # ------------------------------------------------------------------ host

    def run(self):
        self.prepare_hardware()
        self.run_shots()
        self.report()

    def analyse(self):
        """Per-leg delays and statistics from the raw kernel arrays."""
        n = self.n_shots
        stamps = np.asarray(self.stamp_mu, dtype=np.int64).reshape(_N_LEGS, n)
        extras = np.asarray(self.extra_edges,
                            dtype=np.int32).reshape(_N_LEGS, n)
        fires = np.asarray(self.fire_mu,
                           dtype=np.int64).reshape(_N_OUTPUTS, n)

        results = []
        for leg, direction, output, input_name, row, connection in _LEGS:
            stamp = stamps[leg - 1]
            seen = stamp != _NO_EDGE
            delay_mu = np.where(seen, stamp - fires[row], 0).astype(np.int64)
            valid = delay_mu[seen]
            result = {
                "leg": leg,
                "direction": direction,
                "output": output,
                "input": input_name,
                "connection": connection,
                "n_shots": n,
                "n_valid": int(valid.size),
                "edge_seen": seen,
                "delay_mu": delay_mu,
                "extra_edges": extras[leg - 1].copy(),
                "n_extra": int(extras[leg - 1].sum()),
                "mean_mu": float("nan"),
                "std_mu": float("nan"),
                "min_mu": 0,
                "max_mu": 0,
            }
            if valid.size:
                values = valid.astype(np.float64)
                result["mean_mu"] = float(values.mean())
                result["std_mu"] = float(values.std())
                result["min_mu"] = int(valid.min())
                result["max_mu"] = int(valid.max())
            results.append(result)
        return results

    def report(self):
        print("measure_drtio_ttl_skew  n_shots %d" % self.n_shots)
        print("delay = input timestamp - scheduled fire timestamp (ns)")
        print("legs 1-3 compare two node clocks: a negative leg 1-2 delay "
              "is not an error")
        for result in self.analyse():
            self._print_leg(result)
            self._archive_leg(result)

    def _print_leg(self, result):
        head = "leg %d  %s  %s -> %s" % (result["leg"], result["direction"],
                                          result["output"], result["input"])
        if result["n_valid"] == 0:
            print("%s   no edge within %d..%+d ns of the fire time in %d "
                  "shots; expects %s"
                  % (head, -_FIRE_AT_MU, _WINDOW_MU - _FIRE_AT_MU,
                     result["n_shots"], result["connection"]))
            return
        print("%s   mean %.2f ns" % (head, result["mean_mu"]))
        print("       n %d/%d   std %.2f ns   min %d ns   max %d ns   "
              "extra edges %d"
              % (result["n_valid"], result["n_shots"], result["std_mu"],
                 result["min_mu"], result["max_mu"], result["n_extra"]))

    def _archive_leg(self, result):
        prefix = "drtio_ttl_skew_leg%d_" % result["leg"]
        for name in ("delay_mu", "edge_seen", "extra_edges", "mean_mu"):
            self.set_dataset(prefix + name, result[name],
                             broadcast=False, persist=False, archive=True)
