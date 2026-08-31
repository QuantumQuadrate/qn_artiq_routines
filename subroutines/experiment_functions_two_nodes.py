"""Experiments written specifically for true master-satellite two-node use."""


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
