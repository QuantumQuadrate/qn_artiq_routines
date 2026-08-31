"""Native two-node GVS with per-point RTIO-underflow retry."""

from artiq.experiment import EnvExperiment

from GeneralVariableScan_master_satellite_mixin import (
    _CatchUnderflowRetryMixin,
    _GeneralVariableScanMasterSatelliteMixin,
)


class GeneralVariableScan_CatchError_master_satellite_two_nodes(
    _CatchUnderflowRetryMixin,
    _GeneralVariableScanMasterSatelliteMixin,
    EnvExperiment,
):
    """GeneralVariableScan_CatchError_master_satellite_two_nodes

    Run native two-node functions on both master-satellite nodes, retrying
    only a scan point that fails with RTIOUnderflow.
    """

    EXPERIMENT_MODE = "two_nodes"

    def build(self):
        self._build_master_satellite_scan()
        self._build_catch_underflow_arguments()
