"""
This code allows for monitoring the four SPCMs. The counts detected per s can
be viewed with the plot_xyline applets. Any Zotino channels and Urukul channels
that were on before running this code will be left on. For long exposure times
(in ms regime for example) we need to use _counter.gate_rising rather than
ttl.count to avoid overflow.

MASTER-SATELLITE
----------------
The four SPCMs are canonical, master-local detectors on BOTH nodes:

    SPCM_H1 -> Node1 ttl0 / ttl0_counter   (was ttl_SPCM0)
    SPCM_V1 -> Node1 ttl1 / ttl1_counter   (was ttl_SPCM1)
    SPCM_H2 -> Node1 ttl8 / ttl8_counter   (was ttl_SPCM0_OtherNode)
    SPCM_V2 -> Node1 ttl9 / ttl9_counter   (was ttl_SPCM1_OtherNode)

so selected_node does NOT change which detectors are read, and the count-rate
datasets are not node-suffixed either (they are absent from
_single_node_redirect_names), so the same four dataset names are written
whichever node is selected. What selected_node does change is which crate's
DDS/Zotino/TTL hardware base.initialize_hardware() touches, and which node's
ExperimentVariables are loaded -- including t_SPCM_exposure, see prepare().

Base still has to run in single_node mode: the *_rate_dataset name attributes
used below are published by _install_single_node_wiring_metadata, which
returns early in two_nodes mode.
"""

from artiq.experiment import *

import logging
import sys, os
cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")

from utilities.BaseExperiment_master_satellite import (
    BaseExperimentMasterSatellite,
    _DatasetRedirectMixin,
)


class MonitorSPCMinApplet(_DatasetRedirectMixin, EnvExperiment):
    """MonitorSPCMinApplet

    Monitor the four master-local SPCMs and publish their count rates.
    """

    VALID_NODES = ("Node1", "Node2")

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperimentMasterSatellite(experiment=self)
        self.base.build()

        # The detectors are master-local either way; this selects whose
        # hardware is initialized and whose ExperimentVariables are loaded.
        self.setattr_argument(
            "selected_node",
            EnumerationValue(self.VALID_NODES),
            "Node selection",
            tooltip="Node1 = alice, Node2 = bob. All four SPCMs are read "
                    "whichever node is selected; this picks which crate is "
                    "initialized and whose variables are loaded.",
        )

        self.setattr_argument("run_time_minutes", NumberValue(1))
        self.setattr_argument("t_SPCM_exposure", NumberValue(0.05))

        self.setattr_argument(
            "create_applets",
            BooleanValue(True),
            "Applets",
            tooltip="Ask the dashboard to show this node's applets, including "
                    "the five SPCM count-rate plots under Optimization. Untick "
                    "to leave the dashboard exactly as it is. Requires the "
                    "applet dock's CCB policy to be 'Create and enable/disable "
                    "applets'.",
        )

        # The standalone Base archived the GUI arguments here with
        # set_datasets_from_gui_args(). The master-satellite Base has no such
        # method, and ARTIQ already stores the submitted arguments in the HDF5
        # under expid, so nothing is lost by dropping it.
        print("build - done")

    def prepare(self):
        """
        performs initial calculations and sets parameter values before
        running the experiment. also sets data filename for now.

        any conversions from human-readable units to machine units (mu) are done here
        """

        # t_SPCM_exposure is a run-local GUI value that shares its name with
        # the persistent node dataset t_SPCM_exposure_<node>. Capture the
        # submitted value before configure_execution loads that dataset over
        # the same attribute, and re-assert it after base.prepare() -- which
        # runs refresh_compatibility_variables() again -- so the GUI wins for
        # this run. run_time_minutes has no such collision.
        submitted_t_SPCM_exposure = self.t_SPCM_exposure

        node = str(self.selected_node)
        if node not in self.VALID_NODES:
            raise ValueError(
                f"Unsupported selected_node {self.selected_node!r}; expected "
                "'Node1' or 'Node2'."
            )

        self.base.configure_execution("single_node", node)
        # The shared record helpers branch on the alice/bob presentation.
        self.which_node = self.base.NODE_LEGACY_NAMES[node]
        self.base.prepare()

        self.t_SPCM_exposure = submitted_t_SPCM_exposure

        # prepare_laser_stabilizer() is deliberately not called: this is a
        # passive monitor and never touches laser_stabilizer or the
        # per-channel stabilizer_AOM_A* objects.

        self.n_steps = int(60*self.run_time_minutes/self.t_SPCM_exposure+0.5)

        # Seeded on the host rather than inside the kernel: set_dataset is a
        # blocking RPC and would spend RTIO slack before the first event.
        self.set_dataset(self.SPCM0_rate_dataset, [0.0], broadcast=True)
        self.set_dataset(self.SPCM1_rate_dataset, [0.0], broadcast=True)
        self.set_dataset(self.SPCM0_OtherNode_rate_dataset, [0.0], broadcast=True)
        self.set_dataset(self.SPCM1_OtherNode_rate_dataset, [0.0], broadcast=True)
        self.set_dataset(self.AllSPCMs_rate_dataset, [0.0], broadcast=True)

        self._create_node_applets()

        print(self.n_steps)
        print("prepare - done")

    def _create_node_applets(self):
        """Ask the dashboard for this node's applets, if enabled.

        Never fatal: the applets are a convenience, and the CCB is a
        dashboard-side service that artiq_run and the offline compile check
        do not provide meaningfully. Base is always single_node here, so the
        mode guard the GVS mixin needs is unnecessary.
        """
        if not self.create_applets:
            return
        try:
            from applets_master_satellite import create_applets_for

            created = create_applets_for(self, self.base)
        except Exception as error:  # noqa: BLE001 - convenience only
            logging.warning("could not create applets: %s", error)
        else:
            logging.info("requested %d applets for %s", len(created),
                         self.base.which_node)

    @kernel
    def run(self):
        # base.initialize_hardware() owns the core reset in the
        # master-satellite stack, and waits for the satellite when Node2 is
        # the selected node. Nothing is turned off: this is a passive monitor.
        self.base.initialize_hardware(turn_off_dds_channels=False,
                                      turn_off_zotinos=False)

        self.core.break_realtime()
        delay(10 * ms)

        for i in range(self.n_steps):
            with parallel:
                self.SPCM_H1_counter.gate_rising(self.t_SPCM_exposure)
                self.SPCM_V1_counter.gate_rising(self.t_SPCM_exposure)
                self.SPCM_H2_counter.gate_rising(self.t_SPCM_exposure)
                self.SPCM_V2_counter.gate_rising(self.t_SPCM_exposure)

            SPCM_H1_counts = self.SPCM_H1_counter.fetch_count()
            SPCM_V1_counts = self.SPCM_V1_counter.fetch_count()
            SPCM_H2_counts = self.SPCM_H2_counter.fetch_count()
            SPCM_V2_counts = self.SPCM_V2_counter.fetch_count()
            AllSPCMs_counts = (SPCM_H1_counts + SPCM_V1_counts
                               + SPCM_H2_counts + SPCM_V2_counts)

            delay(1 * ms)
            self.append_to_dataset(self.SPCM0_rate_dataset,
                                   SPCM_H1_counts / self.t_SPCM_exposure)
            self.append_to_dataset(self.SPCM1_rate_dataset,
                                   SPCM_V1_counts / self.t_SPCM_exposure)
            self.append_to_dataset(self.SPCM0_OtherNode_rate_dataset,
                                   SPCM_H2_counts / self.t_SPCM_exposure)
            self.append_to_dataset(self.SPCM1_OtherNode_rate_dataset,
                                   SPCM_V2_counts / self.t_SPCM_exposure)
            # The only count-rate applet enabled by default ("All SPCMs count
            # rate", applets_master_satellite.py) subscribes to this one.
            self.append_to_dataset(self.AllSPCMs_rate_dataset,
                                   AllSPCMs_counts / self.t_SPCM_exposure)

        print("Experiment finished.")
