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

THERE IS NO NODE SELECTION, and there is nothing for one to select. The
detectors above are the same four whichever node is "chosen"; the five
count-rate datasets are unsuffixed (they are absent from
_single_node_redirect_names), so the same names are written either way; and
run() below touches no per-node device and no per-node value. Base is
configured single_node/Node1 because Node1 IS the crate the detectors live
on, which makes initializing it both correct and sufficient.

It used to offer Node1/Node2, and the option was worse than useless:
selecting Node2 made base.initialize_hardware() wait for the satellite and
initialize the wrong crate's DDS/Zotino/TTL, for a measurement that reads
neither. The one thing it genuinely changed was which t_SPCM_exposure_<node>
loaded -- and that is 0.01 on both nodes, and the GUI value overrides it
anyway (see prepare()).

Nothing here needs two_nodes mode either. The *_rate_dataset name attributes
come from _install_wiring_metadata, which published them only in single_node
mode until 2026-10-07 and now publishes them in both; single_node/Node1 is
kept because it initializes one crate instead of two and never waits on the
fibre.
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

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperimentMasterSatellite(experiment=self)
        self.base.build()

        self.setattr_argument("run_time_minutes", NumberValue(1))
        self.setattr_argument("t_SPCM_exposure", NumberValue(0.05))

        self.setattr_argument(
            "create_applets",
            BooleanValue(True),
            "Applets",
            tooltip="Ask the dashboard to show the SPCM applets, including "
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

        # NODE1, ALWAYS. There is nothing for a node selection to select:
        # all four SPCMs are master-local on both nodes, the five count-rate
        # datasets are unsuffixed, and run() below touches no per-node device
        # or value. Node1 is simply the crate the detectors physically live
        # on, so initializing it is both correct and sufficient.
        self.base.configure_execution("single_node", "Node1")
        # Published for the shared helpers that branch on the alice/bob
        # presentation; Base only sets it itself on the two-node path.
        self.which_node = self.base.NODE_LEGACY_NAMES["Node1"]
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
            from applets_master_satellite import (
                APPLET_SPECS,
                SPCM_MONITOR_APPLET_SPECS,
                create_applets_for,
            )

            # The SPCM count-rate applets are deliberately NOT in the default
            # set: they would otherwise come up on every scan. This experiment
            # is the one that exists to watch them, so it is the one that asks
            # for them.
            created = create_applets_for(
                self, self.base,
                specs=APPLET_SPECS + SPCM_MONITOR_APPLET_SPECS,
            )
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
