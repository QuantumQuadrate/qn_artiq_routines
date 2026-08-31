"""Single-node compatibility GVS with per-point RTIO-underflow retry."""

from artiq.experiment import EnumerationValue, EnvExperiment

from GeneralVariableScan_master_satellite_mixin import (
    _CatchUnderflowRetryMixin,
    _GeneralVariableScanMasterSatelliteMixin,
)


class GeneralVariableScan_CatchError_master_satellite_single_node(
    _CatchUnderflowRetryMixin,
    _GeneralVariableScanMasterSatelliteMixin,
    EnvExperiment,
):
    """GeneralVariableScan_CatchError_master_satellite_single_node

    Run existing single-node experiment functions on selected Node1 or Node2
    hardware, retrying only a scan point that fails with RTIOUnderflow.
    """

    EXPERIMENT_MODE = "single_node"

    def build(self):
        self.setattr_argument(
            "selected_node", EnumerationValue(self.VALID_NODES)
        )
        self._build_master_satellite_scan()
        self._build_catch_underflow_arguments()
