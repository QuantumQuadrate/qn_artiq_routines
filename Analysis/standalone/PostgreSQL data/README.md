# PostgreSQL data

Monitoring and database tooling: the magnetometer / Zotino / coil-monitor
logging, the Grafana dashboard setup, and how to browse the PostgreSQL
database with Postbird.

## These are pre-master-satellite analysis codes

Everything in this folder — and in the rest of `Analysis/standalone/` — was
written and used **before the master-satellite implementation**. It is kept
working, not retired: any data recorded before the changeover can still be
analysed with these notebooks exactly as it always was.

Data taken *after* the changeover is a different shape. Results are now
written per node, and single-node runs carry the node in the result filename
(`..._Node1_...`, `..._Node2_...`), while master-satellite datasets are
node-suffixed (`AllSPCMs_RO1_Node1`) with a few deliberately kept global
(`n_measurements`, `two_atom_threshold`). Analysis for that lives in
`Analysis/master_satellite/`.

So: old data here, new data there. If a notebook in this folder is ever
pointed at a post-changeover run and comes back empty, the node suffix is the
first thing to check.
