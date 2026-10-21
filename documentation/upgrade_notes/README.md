# Upgrade notes

One file per planned change, named for its **topic** — not dated like
`debug_notes/`, because these are living documents that get edited as the
thinking firms up and as the change is executed in stages.

These are changes that are **understood but deliberately not started**. A note
belongs here when the investigation is finished but the work is deferred, and
the cost of losing the investigation is high — the alias maps, the per-node
hardware differences and the ordering constraints in this codebase take hours
to re-derive and are not visible from any single file.

A good note records, at minimum:

* **what the change is**, concretely enough to execute without re-reading the
  codebase — file paths, line numbers, counts;
* **decisions already taken**, so they are not re-litigated;
* **what makes it safe**, or what makes it dangerous — the specific mechanism
  that means a given edit is or is not a behaviour change;
* **what must be true before starting** — prerequisites, blocking conflicts,
  and anything that needs a human decision first;
* **latent bugs noticed on the way** that are deliberately *not* being fixed,
  so the change stays mechanical and a regression is unambiguous.

The distinction from `debug_notes/`: those are post-mortems of problems that
took real work to diagnose, including the explanations that turned out to be
wrong. These are forward-looking plans for work not yet done.

| Note | Topic | Status |
| --- | --- | --- |
| [D1_optical_pumping_merge_ideas.md](D1_optical_pumping_merge_ideas.md) | Rename `GRIN1and2_dds`→`dds_GRIN1`, `dds_D1_pumping_DP`→`dds_GRIN2`; then give Node1 a real D1 TTL so it can share Node2's OP/EXC code | Partly applied (2026-10-20): alias renames done (GUI names and the hardware check pending); Node1's D1 on ttl15 with an external RF source, Node1 sequences switched; Node1's GRIN1 path and the merge not started |
| [chopped_RO_ideas.md](chopped_RO_ideas.md) | Low-loss fluorescence readout by chopping the FORT out of phase with the readout light; why the DP AOM currently caps the duty cycle at 26% | Idea, partially built |

Both notes above depend on the **same** decision: whether two-node *standalone*
synchronization is retired, which is what frees `ttl12` and `ttl15` as output
TTLs on both nodes. Settle that once and both unblock. Answered 2026-10-19:
retired, cables removed; the D1 note is now partly applied.
