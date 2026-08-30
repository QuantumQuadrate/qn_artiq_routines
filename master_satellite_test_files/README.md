# Master-Satellite Hardware Reference

These files are a previously generated Node1-master / Node2-satellite
ARTIQ configuration from 2026-07-15.

For the current migration work, treat these files as the assumed hardware
reference unless explicitly told otherwise.

Assumptions:

- Node1 = master
- Node2 = satellite
- Satellite DRTIO destination = 1
- The unified device_db is assumed to match the corresponding gateware
  configuration for the purpose of software development.
- Hardware correctness will be verified separately on the experiment system.

Files:

- node1_master_kasli_soc_20260715.json
- node2_satellite_kasli_soc_20260715.json
- device_db_nodes12_20260715.py

Important:

- Do not regenerate gateware unless explicitly requested.
- Do not modify these reference files unless explicitly requested.
- Use these files to determine the low-level master/satellite device mapping
  while implementing the new software architecture.
- Experiment-facing master-satellite names should use the node suffix
  convention, e.g. `dds_FORT_Node1`, `dds_FORT_Node2`.