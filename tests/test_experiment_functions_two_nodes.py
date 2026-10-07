"""Structural guards on subroutines/experiment_functions_two_nodes.py.

AST checks rather than behavioural ones, because the mistakes they catch
COMPILE CLEANLY and run against the wrong hardware, so neither the
hardware-free suite nor the offline compile sweep can see them.

In two-node mode EVERY device names its node -- Base publishes no bare device
aliases, so a bare device name cannot compile. The invariants worth pinning:

  1. A VALUE-CARRYING call must name a node. set(), set_mu(), set_att() and
     set_dac() each carry a frequency, amplitude or voltage, and the two nodes
     disagree about nearly all of them -- f_FORT 245 vs 240 MHz, every coil
     voltage.

  2. A bare name must be one Base actually publishes. The allowed set is
     DERIVED by building a real two-node Base, so it cannot drift: since no
     device aliases are published, any bare device name fails here.

  3. RAW physical names stay out. self.ttl0, self.sampler0 and self.zotino0 all
     RESOLVE in two-node mode and silently point at the master, because
     bind_physical_device calls setattr_device(unified_name) before applying
     the presentation name and Node1's unified names equal the standalone ones.

  4. At most four timeline events may share a timestamp. The SED has 8 FIFOs
     and DISCARDS the overflow with no exception, so this one is about events
     that never happen rather than events that happen wrongly.

  5. The two readouts must switch identical beams, or the retention ratio they
     form is meaningless.

Import-order note: this module imports the two-node functions, which import
artiq names. test_base_experiment_master_satellite installs a minimal stub
lacking TFloat and friends, so load order cannot be relied on; this file
widens the shared stub itself.
"""

import ast
import sys
import types
import unittest
from pathlib import Path


_STUB_MARKER = "_qn_hardware_free_stub"


def _identity_decorator(function=None, **kwargs):
    if function is not None:
        return function
    return lambda decorated: decorated


def _artiq_experiment_stub_or_none():
    """Return the artiq.experiment stub to (re)configure, or None.

    Never shadows a real artiq installation: the ARTIQ repository scan imports
    test files inside a live worker.
    """
    existing_artiq = sys.modules.get("artiq")
    if existing_artiq is not None:
        if getattr(existing_artiq, _STUB_MARKER, False):
            return sys.modules["artiq.experiment"]
        return None
    try:
        import artiq  # noqa: F401
    except ImportError:
        artiq_module = types.ModuleType("artiq")
        setattr(artiq_module, _STUB_MARKER, True)
        experiment_module = types.ModuleType("artiq.experiment")
        artiq_module.experiment = experiment_module
        sys.modules["artiq"] = artiq_module
        sys.modules["artiq.experiment"] = experiment_module
        return experiment_module
    return None


_experiment_stub = _artiq_experiment_stub_or_none()
if _experiment_stub is not None:
    _stub_exports = {
        "EnvExperiment": object,
        "kernel": _identity_decorator,
        "rpc": _identity_decorator,
        "delay": lambda duration: None,
        "delay_mu": lambda duration: None,
        "at_mu": lambda timestamp: None,
        "now_mu": lambda: 0,
        "parallel": None,
        "sequential": None,
        "NumberValue": lambda value, **kwargs: value,
        "BooleanValue": lambda value=False, **kwargs: value,
        "StringValue": lambda value, **kwargs: value,
        "EnumerationValue": lambda values, **kwargs: tuple(values)[0],
        "TBool": bool,
        "TFloat": float,
        "TInt32": int,
        "TInt64": int,
        "TStr": str,
        "MHz": 1e6,
        "kHz": 1e3,
        "ms": 1e-3,
        "us": 1e-6,
        "ns": 1e-9,
        "s": 1.0,
        "V": 1.0,
    }
    for _name, _value in _stub_exports.items():
        if not hasattr(_experiment_stub, _name):
            setattr(_experiment_stub, _name, _value)
    if hasattr(_experiment_stub, "__all__"):
        _experiment_stub.__all__ = sorted(
            set(_experiment_stub.__all__) | set(_stub_exports)
        )


from GeneralVariableScan_master_satellite_mixin import (  # noqa: E402
    build_two_node_function_registry,
)
from utilities.BaseExperiment_master_satellite import (  # noqa: E402
    BaseExperimentMasterSatellite,
)


_MODULE_PATH = Path("subroutines/experiment_functions_two_nodes.py")

#: Calls that carry a per-node value and therefore must name a node.
_VALUE_CARRYING_CALLS = ("set", "set_mu", "set_att", "set_dac")

#: Raw physical device names. Each RESOLVES in two-node mode and points at a
#: specific crate while reading as node-neutral.
_RAW_DEVICE_PREFIXES = ("ttl", "sampler", "zotino", "urukul", "spi_")

#: Non-device bare names the file legitimately reads: run-global state, host
#: result scalars and buffers, and the dataset helpers.
_ALLOWED_NON_DEVICE = {
    "core", "core_dma", "scheduler",
    "measurement", "measurements_progress", "n_measurements",
    "two_atom_threshold", "two_atom_threshold_for_loading",
    "t_Node2_rtio_offset_mu", "t_Node2_excitation_delay_mu",
    "max_atom_check_tries_two_node", "max_loading_rounds_two_node",
    "n_alternating_RO_windows_per_node",
    "t_alternating_RO_window", "t_alternating_RO_pad",
    "atom_loading_time", "atom_loading_time_list",
    "n_atom_loaded_per_iteration", "n_feedback_per_iteration",
    "in_health_check",
    "AllSPCMs_RO1", "AllSPCMs_RO2",
    "AllSPCMs_RO1_list", "AllSPCMs_RO2_list",
    "AllSPCMs_alternating_RO_alice", "AllSPCMs_alternating_RO_bob",
    "print_async", "append_to_dataset", "set_dataset", "get_dataset",
    "write_results",
}


def _self_attributes(node):
    """Every attribute read off `self`, with its line number."""
    return [
        (item.attr, item.lineno)
        for item in ast.walk(node)
        if isinstance(item, ast.Attribute)
        and isinstance(item.value, ast.Name)
        and item.value.id == "self"
    ]


class ExperimentFunctionsTwoNodesStructureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = _MODULE_PATH.read_text(encoding="utf-8")
        cls.tree = ast.parse(cls.source)
        cls.self_attributes = _self_attributes(cls.tree)

        # DERIVED, not hand-listed: actually build a two-node Base against the
        # hardware-free fake and collect every attribute it ends up with. That
        # is by construction the set of names a two-node kernel can reach, so
        # the guard cannot drift out of date as Base grows -- and if Base ever
        # stops publishing one of them, this test fails for the right reason
        # rather than silently widening.
        import test_base_experiment_master_satellite as base_tests

        experiment = base_tests.FakeExperiment()
        base = BaseExperimentMasterSatellite(experiment, "two_nodes")
        base.build()
        base.prepare()
        base.initialize_result_state()

        cls.published_bare = {
            name for name in vars(experiment)
            if not name.endswith(("_Node1", "_Node2"))
        } | _ALLOWED_NON_DEVICE

    def test_module_touches_self_attributes_at_all(self):
        """Guard the guard: if the walk finds nothing, the rest proves nothing."""
        self.assertGreater(len(self.self_attributes), 100)

    def test_value_carrying_calls_always_name_a_node(self):
        """set()/set_mu()/set_att()/set_dac() must target a _NodeX device.

        This is THE invariant of the broadcast design. A bare .set() would
        carry one frequency to two crates that disagree about it. The
        broadcast objects have no set() so the compiler already refuses, but
        set_dac on the Zotinos is not covered by that -- the Zotinos are not
        broadcast -- and this catches both in one place.
        """
        offenders = []
        for node in ast.walk(self.tree):
            if not isinstance(node, ast.Call):
                continue
            func = node.func
            if not isinstance(func, ast.Attribute):
                continue
            if func.attr not in _VALUE_CARRYING_CALLS:
                continue
            target = func.value
            # self.<device>.set(...)
            if not (isinstance(target, ast.Attribute)
                    and isinstance(target.value, ast.Name)
                    and target.value.id == "self"):
                continue
            if not target.attr.endswith(("_Node1", "_Node2")):
                offenders.append(
                    f"{_MODULE_PATH}:{func.lineno}: "
                    f"self.{target.attr}.{func.attr}(...)"
                )
        self.assertEqual(
            offenders, [],
            "These carry a per-node value but address a bare (broadcast) "
            "device. The two nodes disagree about nearly every frequency, "
            "amplitude and coil voltage, so each must name its node:\n  "
            + "\n  ".join(offenders),
        )

    def test_bare_names_are_ones_base_actually_publishes(self):
        """A bare name must be a broadcast alias, a shared detector, or an
        agreeing scalar -- not a per-node name someone forgot to suffix."""
        offenders = []
        for name, line in self.self_attributes:
            if name.endswith(("_Node1", "_Node2")) or name.startswith("_"):
                continue
            if name in self.published_bare:
                continue
            offenders.append(f"{_MODULE_PATH}:{line}: self.{name}")
        self.assertEqual(
            offenders, [],
            "Not published bare in two-node mode. Either suffix it with the "
            "node, or add it to TWO_NODE_AGREEING_SCALARS if both nodes "
            "genuinely share one value:\n  " + "\n  ".join(offenders),
        )

    def test_no_raw_physical_device_names(self):
        """self.ttl0 / self.sampler0 / self.zotino0 resolve -- to the MASTER."""
        offenders = []
        for name, line in self.self_attributes:
            if name.endswith(("_Node1", "_Node2")):
                continue
            for prefix in _RAW_DEVICE_PREFIXES:
                if not name.startswith(prefix):
                    continue
                tail = name[len(prefix):]
                if tail and (tail[0].isdigit() or prefix == "spi_"):
                    offenders.append(f"{_MODULE_PATH}:{line}: self.{name}")
        self.assertEqual(
            offenders, [],
            "Raw physical device names silently address one crate. Use the "
            "broadcast alias, or a _NodeX name:\n  " + "\n  ".join(offenders),
        )

    def test_registry_contains_exactly_the_intended_experiments(self):
        """The registry matches on the substring "experiment", so a helper
        that picked it up would appear in the dashboard as runnable."""
        registry = build_two_node_function_registry()
        self.assertEqual(
            set(registry),
            {
                "master_satellite_namespace_sanity_experiment",
                "Two_nodes_atom_loading_experiment",
                "Two_nodes_alternating_shot_experiment",
            },
        )
        module_functions = [
            node.name for node in self.tree.body
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        ]
        for name in module_functions:
            if name in registry or name.startswith("_"):
                continue
            self.assertNotIn(
                "experiment", name,
                f"{name} would be picked up by the registry; rename it or "
                f"make it private.",
            )

    def test_every_gate_is_matched_by_exactly_one_fetch_per_counter(self):
        """Each gate_rising pushes one total; each fetch_count pops one.

        Fetching twice for one gate reads a total that is not there; fetching
        too few leaves the FIFO filling. Both compile cleanly.
        """
        for node in self.tree.body:
            if not isinstance(node, ast.FunctionDef):
                continue
            gates = fetches = 0
            for inner in ast.walk(node):
                if not isinstance(inner, ast.Call):
                    continue
                func = inner.func
                if not isinstance(func, ast.Attribute):
                    continue
                target = func.value
                if not (isinstance(target, ast.Attribute)
                        and target.attr.endswith("_counter")):
                    continue
                if func.attr.startswith("gate_rising"):
                    gates += 1
                elif func.attr == "fetch_count":
                    fetches += 1
            if gates == 0 and fetches == 0:
                continue
            self.assertEqual(
                gates, fetches,
                f"{node.name} opens {gates} counter gate(s) but performs "
                f"{fetches} fetch_count call(s). Each gate must be read "
                f"exactly once -- an edge counter is drained by reading it.",
            )

    # artiq/gateware/rtio/sed/core.py: lane_count=8.
    SED_LANE_COUNT = 8
    # Lab rule: group timeline events in fours, so there is 2x headroom on the
    # 8 FIFOs rather than the zero margin an 8-event group would have.
    MAX_EVENTS_PER_TIMESTAMP = 4

    def test_no_run_of_switch_events_shares_a_timestamp(self):
        """Consecutive switch calls with no delay between them are DROPPED.

        Every .sw.on()/.off() is one rtio_output. From the SED design note,
        artiq/gateware/rtio/sed/__init__.py:

            "When an event is submitted, it is written into the current FIFO
             if its timestamp is strictly increasing. Otherwise, the current
             FIFO number is incremented by one (and wraps around, if the
             current FIFO was the last) and the event is written there, unless
             that FIFO already contains an event with a greater timestamp. In
             that case, an asynchronous error is reported."

            "The maximum number of simultaneous events (on different
             channels), and the maximum number of active timeline 'rewinds',
             are equal to the number of FIFOs."

        So the budget is 8, and SIMULTANEOUS EVENTS AND REWINDS SHARE IT. Past
        8, do_write stays 0 and the event is silently discarded: no exception,
        just an "RTIO sequence error involving channel ..." line in the core
        log and a switch that never moved.

        This guard counts straight-line runs of switch events only. Rewinds --
        a backward at_mu, or each statement after the first inside `with
        parallel:` -- draw on the same 8. The file's largest `with parallel:`
        gates four counters, so there is headroom, but a future change that
        adds rewinds next to a long run can still exceed the budget without
        tripping this test.

        Observed on hardware 2026-10-07 (RID 38676): the loader turned on both
        repumps, both cooling DPs and all twelve fiber AOMs at one timestamp --
        16 events. Events 9-16 were dropped, so A3/A4/A5/A6 never switched on
        EITHER node and the MOT ran on two of six beams. The log named
        0x01002f = ttl_urukul4_sw1 = dds_AOM_A6_Node2, the 16th and last event.

        Two-node mode is where this bites: the legacy two-node code ran one
        kernel per crate, so each issued only 8. Merging them onto one timeline
        doubled every switching point.

        The rule, from the lab owner: at least 5 us between consecutive switch
        events. That also clears the separate slack-erosion failure, since the
        cursor then advances 5 us per event against ~0.7 us of submission cost.
        """
        def is_event(stmt):
            """Anything that puts an event on the timeline without advancing it."""
            if isinstance(stmt, ast.Expr) and isinstance(stmt.value, ast.Call):
                func = stmt.value.func
                if isinstance(func, ast.Attribute):
                    if func.attr in ("on", "off") and not stmt.value.args:
                        return True
                    # dds.set / zotino.set_dac / dds.set_att count too: they are
                    # timeline events like any other.
                    if func.attr in ("set", "set_dac", "set_att"):
                        return True
            return False

        def advances_cursor(stmt):
            # NOT ast.With: `with parallel:` restarts every branch at the
            # block's entry time, so each branch's FIRST event shares one
            # timestamp with every other branch's first event. Treating the
            # block as opaque would hide exactly the collision this guard
            # exists to find. Handled explicitly in walk() instead.
            if False:
                return True
            if isinstance(stmt, ast.Expr) and isinstance(stmt.value, ast.Call):
                func = stmt.value.func
                name = (
                    func.id if isinstance(func, ast.Name)
                    else getattr(func, "attr", "")
                )
                # at_mu COUNTS as advancing. That is an assumption, and it
                # holds because every at_mu left in this file targets a
                # computed FORWARD offset (the alternating readout's window
                # grid). What would break it is the no-op pattern -- capture
                # `x = now_mu()` while the cursor still sits on an event's
                # timestamp, then at_mu(x) -- which separates nothing. The
                # per-node stages that used to do exactly that now use
                # `with parallel` instead, and
                # test_no_now_mu_anchor_captures_a_live_timestamp keeps it so.
                return name in ("delay", "delay_mu", "at_mu", "break_realtime")
            return False

        def first_statements(stmts):
            """Statements that can execute first -- an unguarded `if` may be
            skipped, so it falls through to whatever follows it."""
            reachable = []
            for stmt in stmts:
                if isinstance(stmt, ast.If):
                    reachable += first_statements(stmt.body)
                    if stmt.orelse:
                        reachable += first_statements(stmt.orelse)
                        break
                    continue
                reachable.append(stmt)
                break
            return reachable

        def first_statements_that_are_events(stmts):
            """Branch-entry statements that put an event on the timeline.

            In a `with parallel:` block each top-level statement is its own
            branch and every branch restarts at the block's entry time, so
            these all land on ONE timestamp.
            """
            return [
                stmt
                for branch in stmts
                for stmt in first_statements([branch])
                if is_event(stmt)
            ]

        offenders = []

        def walk(stmts, counts, function_name):
            """Propagate the running event count along every path."""
            for stmt in stmts:
                if advances_cursor(stmt):
                    counts = {0}
                    continue
                if is_event(stmt):
                    counts = {c + 1 for c in counts}
                    worst = max(counts)
                    if worst > self.MAX_EVENTS_PER_TIMESTAMP:
                        offenders.append((function_name, stmt.lineno, worst))
                    continue
                if isinstance(stmt, ast.If):
                    taken = walk(stmt.body, set(counts), function_name)
                    if stmt.orelse:
                        counts = taken | walk(
                            stmt.orelse, set(counts), function_name
                        )
                    else:
                        counts = taken | counts
                    continue
                if isinstance(stmt, ast.With):
                    # Every branch starts at the block's entry timestamp, so
                    # their first events collide with each other and with
                    # whatever already shares that timestamp.
                    entry = max(counts)
                    starters = len(first_statements_that_are_events(stmt.body))
                    if entry + starters > self.MAX_EVENTS_PER_TIMESTAMP:
                        offenders.append(
                            (function_name, stmt.lineno, entry + starters)
                        )
                    outgoing = set()
                    for branch in stmt.body:
                        outgoing |= walk([branch], set(counts), function_name)
                    counts = outgoing or {0}
                    continue
                if isinstance(stmt, (ast.For, ast.While)):
                    exits = walk(stmt.body, set(counts), function_name)
                    # a second pass with the loop's own exit counts, so a
                    # group that straddles the back-edge is seen
                    walk(stmt.body, set(exits), function_name)
                    counts = exits | counts
            return counts

        for node in self.tree.body:
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                walk(node.body, {0}, node.name)

        self.assertFalse(
            offenders,
            "more than "
            + str(self.MAX_EVENTS_PER_TIMESTAMP)
            + " timeline events share one timestamp (function, line, count): "
            + repr(offenders)
            + f". The gateware DISCARDS silently past {self.SED_LANE_COUNT}. "
            f"Put delay(5 * us) after every "
            f"{self.MAX_EVENTS_PER_TIMESTAMP}.",
        )

    def test_the_two_readouts_set_up_identical_optics(self):
        """RO1 and RO2 must leave the beams in the same state.

        Retention is RO2 compared against RO1, so any beam that is on for one
        shot and off for the other does not merely add noise -- it makes the
        ratio meaningless.

        Found on hardware 2026-10-07 from "readout1 and readout2 have different
        fluorescence distribution": PGC_and_RO_with_on_chip_beams is True on
        both nodes, second_shot had an `else` turning A5/A6 OFF, and first_shot
        had no `else` at all -- so A5/A6 kept the ON state the loader left them
        in. The first readout ran on six beams and the second on four. The same
        asymmetry is still present in the single-node original.

        Compares the per-node-suffix-stripped action shapes, so it checks the
        STRUCTURE (which device, which branch condition) rather than literal
        text.
        """
        functions = {
            node.name: node for node in self.tree.body
            if isinstance(node, ast.FunctionDef)
        }

        def optical_actions(function):
            actions = []

            def walk(stmts, conditions):
                for stmt in stmts:
                    if isinstance(stmt, ast.If):
                        test = ast.unparse(stmt.test)
                        walk(stmt.body, conditions + [test])
                        if stmt.orelse:
                            walk(stmt.orelse, conditions + ["NOT " + test])
                        continue
                    if isinstance(stmt, ast.With):
                        walk(stmt.body, conditions)
                        continue
                    if (isinstance(stmt, ast.Expr)
                            and isinstance(stmt.value, ast.Call)
                            and isinstance(stmt.value.func, ast.Attribute)):
                        func = stmt.value.func
                        if func.attr in ("on", "off"):
                            target = ast.unparse(func.value)
                            actions.append((
                                (target + "." + func.attr).replace(
                                    "_Node1", "_N").replace("_Node2", "_N"),
                                tuple(
                                    c.replace("_Node1", "_N").replace(
                                        "_Node2", "_N")
                                    for c in conditions
                                ),
                            ))

            walk(function.body, [])
            counts = {}
            for action in actions:
                counts[action] = counts.get(action, 0) + 1
            return counts

        first = optical_actions(functions["first_shot"])
        second = optical_actions(functions["second_shot"])

        only_first = {k: v for k, v in first.items() if second.get(k) != v}
        only_second = {k: v for k, v in second.items() if first.get(k) != v}

        self.assertEqual(
            (only_first, only_second), ({}, {}),
            "first_shot and second_shot switch different beams, so RO1 and RO2 "
            "are not comparable and the retention ratio is meaningless. "
            f"only in first_shot: {only_first}; only in second_shot: "
            f"{only_second}",
        )

    def test_no_now_mu_anchor_captures_a_live_timestamp(self):
        """`x = now_mu()` must not sit on an event's timestamp.

        The pattern that breaks the same-timestamp guard is: emit an event,
        capture `x = now_mu()` while the cursor is still ON that event's
        timestamp, then at_mu(x) -- which separates nothing, so every event
        placed at that anchor piles onto the original. It reads like clean
        anchored code and the guard's "at_mu advances" assumption hides it.

        Found by an audit on 2026-10-07 as a 3-event pile-up in the loader's
        PGC block. The fix was structural rather than a pad: `with parallel`
        gives a common start and a resync past the longest branch for free, so
        the per-node stages need no anchor at all.

        Allowed: a capture preceded by something that ADVANCES the cursor, so
        the anchor lands on fresh time.
        """
        def is_event(stmt):
            return (
                isinstance(stmt, ast.Expr)
                and isinstance(stmt.value, ast.Call)
                and isinstance(stmt.value.func, ast.Attribute)
                and stmt.value.func.attr in ("on", "off", "set", "set_dac",
                                             "set_att")
            )

        def is_now_mu_capture(stmt):
            return (
                isinstance(stmt, ast.Assign)
                and isinstance(stmt.value, ast.Call)
                and isinstance(stmt.value.func, ast.Name)
                and stmt.value.func.id == "now_mu"
            )

        offenders = []

        def scan(stmts, function_name):
            previous = None
            for stmt in stmts:
                if is_now_mu_capture(stmt) and previous is not None:
                    if is_event(previous):
                        offenders.append((function_name, stmt.lineno))
                for attribute in ("body", "orelse", "finalbody"):
                    inner = getattr(stmt, attribute, None)
                    if isinstance(inner, list) and inner:
                        scan(inner, function_name)
                previous = stmt

        for node in self.tree.body:
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                scan(node.body, node.name)

        self.assertFalse(
            offenders,
            "now_mu() captured while the cursor still sits on an event's "
            f"timestamp (function, line): {offenders}. Every at_mu to that "
            f"anchor piles onto that event. Put a delay before the capture, "
            f"or use `with parallel` and drop the anchor.",
        )

    def test_every_amplitude_index_read_is_an_index_fed_back(self):
        """A setpoint index the sequence READS must be one feedback WRITES.

        FeedbackChannel.amplitudes is np.zeros(len(set_points)) with only [0]
        seeded (aom_feedback.py:106), and run(setpoint_index=N) writes ONLY
        amplitudes[N]. An index that is never fed back therefore stays 0.0 --
        and 0.0 read back as a dds amplitude turns the FORT OFF rather than
        raising, so nothing anywhere complains.

        Found on 2026-10-07: the interleaved two-node feedback paired one node
        with one index (Node1 at 2, Node2 at 1) instead of running each node at
        both, so stabilizer_FORT_Node1.amplitudes[1] stayed 0.0 while
        first_shot, second_shot and two_node_alternating_shot all read exactly
        that index. Node1's FORT would have been commanded to zero amplitude at
        every readout.
        """
        # indices written: stabilizer_<x>_<Node>.run(setpoint_index=N), with a
        # bare run() meaning index 0
        written = {}
        # indices read: self.stabilizer_<x>_<Node>.amplitudes[N]
        read = {}

        for inner in ast.walk(self.tree):
            if isinstance(inner, ast.Call):
                func = inner.func
                if (isinstance(func, ast.Attribute) and func.attr == "run"
                        and isinstance(func.value, ast.Attribute)
                        and func.value.attr.startswith("stabilizer_")):
                    channel = func.value.attr
                    index = 0
                    for keyword in inner.keywords:
                        if keyword.arg == "setpoint_index":
                            self.assertIsInstance(
                                keyword.value, ast.Constant,
                                f"{channel}.run() setpoint_index must be a "
                                f"literal for this guard to reason about it",
                            )
                            index = keyword.value.value
                    written.setdefault(channel, set()).add(index)
            if isinstance(inner, ast.Subscript):
                value = inner.value
                if (isinstance(value, ast.Attribute)
                        and value.attr == "amplitudes"
                        and isinstance(value.value, ast.Attribute)
                        and value.value.attr.startswith("stabilizer_")):
                    if isinstance(inner.slice, ast.Constant):
                        read.setdefault(value.value.attr, set()).add(
                            inner.slice.value
                        )

        self.assertTrue(read, "no amplitudes[...] reads found; guard is inert")

        for channel, indices in sorted(read.items()):
            # laser_stabilizer_<Node>.run() also feeds back the FORT channel at
            # index 0, so index 0 is always available.
            fed_back = written.get(channel, set()) | {0}
            missing = sorted(indices - fed_back)
            self.assertFalse(
                missing,
                f"{channel}.amplitudes{missing} is read by the sequence but "
                f"never written: run() is only called at setpoint_index "
                f"{sorted(written.get(channel, set()))}. Those entries are "
                f"still np.zeros, so the dds would be set to amplitude 0.0.",
            )

    def test_gating_always_uses_all_four_canonical_detectors(self):
        """All four SPCMs see both nodes, so a gate must open all four.

        Forgetting one reads a quarter low and still runs.
        """
        expected = {
            "ttl_SPCM0_counter", "ttl_SPCM1_counter",
            "ttl_SPCM0_OtherNode_counter", "ttl_SPCM1_OtherNode_counter",
        }
        for node in self.tree.body:
            if not isinstance(node, ast.FunctionDef):
                continue
            gated = {
                func.value.attr
                for func in (
                    item.func for item in ast.walk(node)
                    if isinstance(item, ast.Call)
                )
                if isinstance(func, ast.Attribute)
                and func.attr.startswith("gate_rising")
                and isinstance(func.value, ast.Attribute)
            }
            if not gated:
                continue
            self.assertEqual(
                gated, expected,
                f"{node.name} gates {sorted(gated)}; it must gate all four "
                f"canonical detectors.",
            )

    def test_sanity_attributes_are_all_reachable_concepts(self):
        """The sanity list must name only things two-node mode really has."""
        from subroutines import experiment_functions_two_nodes as module
        for name in module.MASTER_SATELLITE_SANITY_ATTRIBUTES:
            if name in self.published_bare:
                continue
            self.assertTrue(
                name.endswith(("_Node1", "_Node2")),
                f"{name} is neither node-suffixed nor published bare, so the "
                f"sanity check would look for something that cannot exist.",
            )


if __name__ == "__main__":
    unittest.main()
