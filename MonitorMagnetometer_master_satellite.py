"""
Monitor one magnetometer, with one node's hardware running.

Ported from tests/A6-MonitorMagnetometer.py, which stays standalone and
unchanged. The class name still contains "MonitorMagnetometer" so that
Analysis/master_satellite/2026_Monitor_Magnetometer_Analysis.ipynb, whose
filter is name_filters=["MonitorMagnetometer"], picks these runs up without
being touched.

MASTER-SATELLITE
----------------
Two independent selections, because the node whose hardware is running and
the magnetometer being read do not have to be the same one:

  selected_node       Base runs in single_node mode for this node, so
                      zotino0, coil_channels, dds_FORT and AZ_bottom_volts_MOT
                      and friends all resolve to its crate. The OTHER node is
                      driven dark by base.force_other_node_off().

  magnetometer_node   Which node's magnetometer is sampled. Normally the same
                      node; setting it to the other one measures that node's
                      field while its own beams and coils are off, which is
                      the cross-talk measurement (how much of this node's
                      coil current the other trap sees).

SAMPLER PORTS
-------------
Each node has its own magnetometer and each sits on that node's own sampler2,
on channels 1/2/3 for the sensor's X/Y/Z axes:

    node    node-local   physical    RTIO channels        crate
    Node1   sampler2     sampler2    spi 30,    cnv 32    master     (dest 0)
    Node2   sampler2     sampler5    spi 65566, cnv 65568 satellite  (dest 1)

The physical mapping is utilities/config/master_satellite/device_aliases.json;
the channel indices are BaseExperimentMasterSatellite's
SINGLE_NODE_WIRING_METADATA, which declares Magnetometer_X_ch/Y_ch/Z_ch as
1/2/3 for BOTH nodes. prepare() reads them out of the MAGNETOMETER's node
rather than using the bare attributes single_node mode publishes for the
selected node, so the two selections stay independent if that wiring ever
diverges. (The standalone docstring says "Sampler2 Ch2,3,4", which is the same
three channels counted from one.)

Gain: _initialize_sampler_group sets explicit gains only on Node1's sampler0
(x8) and sampler1 (x1). sampler2 is left at what init() gives it, gain 0, so
on both nodes the magnetometer is read on the plain +/-10 V range.

WHY THE SAMPLER IS RE-INITIALIZED HERE
--------------------------------------
single_node mode only removes the other node from the INITIALIZATION
lifecycle: _deactivate_node_hardware_groups rebinds _samplers_<other node> to
[], so _initialize_sampler_group never sees that node's sampler2 and never
calls init() on it, even though build() bound the device. Reading it without
init() would sample whatever the SPI core was last left in. So run() inits the
selected magnetometer's sampler explicitly. Unconditional, because re-initing
the selected node's own sampler2 is a plain SPI reconfiguration.

UNITS
-----
V_to_mG, default 1.0, i.e. raw volts by default. The standalone branched on
which_node and the two branches disagreed: alice stored milligauss (x350) and
bob stored volts with "conversion done in analysis code". Now that either
magnetometer can be selected, a node-dependent unit would make the two files
incomparable, so the factor is one explicit argument instead of a hidden
branch, and the default matches what the analysis notebook expects -- it
applies its own V_to_mG = 625 to the stored numbers.

The axis permutation IS kept exactly as the standalone had it for both nodes:
stored Magnetometer_X is the sensor's Y channel and stored Magnetometer_Y is
the sensor's X channel, because the sensor's X axis is the coils' Y axis.
Note this disagrees with measure_Magnetometer() in experiment_functions,
which maps Zero_X <- Z, Zero_Y <- X, Zero_Z <- Y. That inconsistency predates
this file and is left alone here; this experiment's own convention is the one
the magnetometer notebook was written against.

DATASETS
--------
Magnetometer_X/Y/Z, unsuffixed. They are not in Base's
_single_node_redirect_names, so _DatasetRedirectMixin passes them through and
both nodes write the same three names -- which is right, since one run reads
one magnetometer. It also means any magnetometer applet already configured by
hand on the dashboard keeps working unchanged.

Seeded in prepare() rather than in the kernel as the standalone does it:
append_to_dataset can only mutate a dataset that already exists in THIS
worker's dicts, because it resolves through DatasetManager._get_mutation_target
which does not fall through to the master's database the way get() does.
"""

from artiq.experiment import *
import numpy as np

import sys, os
### get the current working directory
current_working_directory = os.getcwd()
cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")

from utilities.BaseExperiment_master_satellite import (
    BaseExperimentMasterSatellite,
    _DatasetRedirectMixin,
)


class MonitorMagnetometer_master_satellite(_DatasetRedirectMixin, EnvExperiment):
    """MonitorMagnetometer_master_satellite

    Monitor either node's magnetometer with one node's hardware running.
    """

    VALID_NODES = ("Node1", "Node2")

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperimentMasterSatellite(experiment=self)
        self.base.build()

        self.setattr_argument(
            "selected_node",
            EnumerationValue(self.VALID_NODES),
            "Node selection",
            tooltip="Node1 = alice, Node2 = bob. The node whose coils, FORT "
                    "and beams are used; the other node is driven dark.",
        )
        self.setattr_argument(
            "magnetometer_node",
            EnumerationValue(self.VALID_NODES),
            "Node selection",
            tooltip="Which node's magnetometer to sample. Each node has its "
                    "own, on that node's sampler2 channels 1/2/3. Set it to "
                    "the node that is NOT selected to measure cross-talk: "
                    "that node's field with its own beams and coils off.",
        )

        self.setattr_argument("n_measurements", NumberValue(100, type='int', scale=1, ndecimals=0, step=1))
        self.setattr_argument("n_average", NumberValue(10, type='int', scale=1, ndecimals=0, step=1))
        self.setattr_argument("t_step_ms", NumberValue(10, type='int', scale=1, ndecimals=0, step=1))
        self.setattr_argument("Coils_settings", EnumerationValue(['MOT', 'Optical_pumping', 'Zero_volts', 'None']))

        self.setattr_argument(
            "V_to_mG",
            NumberValue(1.0, type='float', ndecimals=1, step=1.0),
            tooltip="Factor applied before storing. 1.0 stores raw volts, "
                    "which is what the analysis notebook expects (it applies "
                    "its own V_to_mG = 625). Set 350 or 625 to store "
                    "milligauss directly.",
        )

        self.setattr_argument("enable_individual_coil_scan_mode", BooleanValue(default=False), "individual_coil_scan_mode")
        self.setattr_argument("coil_test_voltage", NumberValue(10, type='int', scale=1, ndecimals=0, step=1))

        ### base.set_datasets_from_gui_args() has no master-satellite
        ### equivalent; the one shadowing argument is re-asserted in prepare().
        print("build - done")

    def prepare(self):
        """
        performs initial calculations and sets parameter values before
        running the experiment. also sets data filename for now.

        any conversions from human-readable units to machine units (mu) are done here
        """
        ### n_measurements is a run-local GUI value that shares its name with a
        ### global persistent variable. Capture it before configure_execution
        ### loads that dataset over the same attribute, and re-assert it after,
        ### so the submitted value wins for this run.
        submitted_n_measurements = self.n_measurements

        node = str(self.selected_node)
        if node not in self.VALID_NODES:
            raise ValueError(
                f"Unsupported selected_node {self.selected_node!r}; expected "
                "'Node1' or 'Node2'."
            )
        magnetometer_node = str(self.magnetometer_node)
        if magnetometer_node not in self.VALID_NODES:
            raise ValueError(
                f"Unsupported magnetometer_node {self.magnetometer_node!r}; "
                "expected 'Node1' or 'Node2'."
            )

        self.base.configure_execution("single_node", node)
        ### The reused single-node code and the shared config paths read
        ### self.which_node expecting the 'alice'/'bob' presentation.
        self.which_node = self.base.NODE_LEGACY_NAMES[node]
        self.base.prepare()

        self.n_measurements = submitted_n_measurements

        self._selected_node = node
        self._magnetometer_node = magnetometer_node

        ### build() bound BOTH nodes' samplers under their suffixed names --
        ### configure_execution ran after base.build(), so _presentation_name
        ### was still suffixing -- and _publish_single_node_physical_
        ### presentations only ADDS the bare aliases for the selected node. So
        ### the other node's sampler2 is still reachable here by its suffixed
        ### name; it is just not initialized, which run() handles.
        sampler_attribute = f"sampler2_{magnetometer_node}"
        try:
            self.magnetometer_sampler = getattr(self, sampler_attribute)
        except AttributeError as error:
            raise RuntimeError(
                f"{magnetometer_node}'s magnetometer sampler "
                f"{sampler_attribute!r} is not bound. Both nodes' samplers are "
                "bound in base.build(); this means execution was configured "
                "before it, which binds only the selected node."
            ) from error

        wiring = self.base.SINGLE_NODE_WIRING_METADATA[magnetometer_node]
        self.magnetometer_x_ch = wiring["Magnetometer_X_ch"]
        self.magnetometer_y_ch = wiring["Magnetometer_Y_ch"]
        self.magnetometer_z_ch = wiring["Magnetometer_Z_ch"]

        ### Seeded here, not in the kernel: append_to_dataset can only mutate
        ### a dataset that already exists in this worker's dicts.
        self.set_dataset("Magnetometer_X", [0.0], broadcast=True)
        self.set_dataset("Magnetometer_Y", [0.0], broadcast=True)
        self.set_dataset("Magnetometer_Z", [0.0], broadcast=True)

        if magnetometer_node != node:
            print(f"NOTE: running {node} but sampling {magnetometer_node}'s "
                  f"magnetometer, on {sampler_attribute}. "
                  f"{magnetometer_node}'s own beams and coils are off.")

        print("prepare - done")

    @kernel
    def run(self):
        ### base.initialize_hardware() owns the core reset in the
        ### master-satellite stack, and waits for the satellite when Node2 is
        ### the selected node.
        self.base.initialize_hardware()

        ### Drive the idle node dark, as the other single-node master-satellite
        ### experiments do. single_node mode only REMOVES the other node from
        ### the initialization lifecycle and emits no hardware operation, so
        ### without this its beams stay on and its coils stay energised from
        ### whatever the previous run left -- which would also put a field on
        ### whichever magnetometer is being read.
        self.base.force_other_node_off()

        self._initialize_magnetometer_sampler()

        self.expt()
        print("*************   Experiment finished   *************")

    @kernel
    def _initialize_magnetometer_sampler(self):
        """init() the sampler the selected magnetometer sits on.

        Needed because _initialize_sampler_group only ever sees the selected
        node's samplers; see WHY THE SAMPLER IS RE-INITIALIZED HERE above.
        break_realtime first: force_other_node_off has just spent real time on
        three CPLD inits and a Zotino init, over remote SPI when the other node
        is the satellite.
        """
        self.core.break_realtime()
        self.magnetometer_sampler.init()
        delay(1 * ms)

    @kernel
    def expt(self):
        """
        The experiment loop.

        :return:
        """
        if not self.enable_individual_coil_scan_mode:
            self.core.break_realtime()

            self.zotino0.set_dac(
                [0.0, 0.0, 0.0, 0.0],
                channels=self.coil_channels)

            delay(1 * ms)
            self.dds_FORT.sw.on()

            if self.Coils_settings == "MOT":
                ### Set the coils to MOT loading setting
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
                    channels=self.coil_channels)
                delay(1 * ms)

            if self.Coils_settings == "Optical_pumping":
                ### Set the coils to OP values
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
                    channels=self.coil_channels)
                delay(1 * ms)

            if self.Coils_settings == "Zero_volts":
                ### Turn off all the coils
                self.zotino0.set_dac(
                    [0.0, 0.0, 0.0, 0.0],
                    channels=self.coil_channels)
                delay(1 * ms)

            if self.Coils_settings == "None":
                ### Does not change the coils values.
                pass
                delay(1 * ms)

            for n_measurement in range(self.n_measurements):

                if n_measurement % max(1, self.n_measurements // 10) == 0:
                    self.print_async("progress (%): ", (n_measurement / self.n_measurements) * 100)

                measurement_buf = np.array([0.0] * 8)
                MagnetometerX = 0.0
                MagnetometerY = 0.0
                MagnetometerZ = 0.0

                for i in range(self.n_average):
                    self.magnetometer_sampler.sample(measurement_buf)
                    MagnetometerX += measurement_buf[self.magnetometer_x_ch]
                    MagnetometerY += measurement_buf[self.magnetometer_y_ch]
                    MagnetometerZ += measurement_buf[self.magnetometer_z_ch]
                    delay(10 * us)
                MagnetometerX /= self.n_average
                MagnetometerY /= self.n_average
                MagnetometerZ /= self.n_average

                ### V_to_mG defaults to 1.0, i.e. raw volts; see UNITS above.
                ### The sensor's X axis is the coils' Y axis, and vice versa,
                ### so the first two are deliberately crossed.
                self.append_to_dataset("Magnetometer_X", MagnetometerY * self.V_to_mG)
                self.append_to_dataset("Magnetometer_Y", MagnetometerX * self.V_to_mG)
                self.append_to_dataset("Magnetometer_Z", MagnetometerZ * self.V_to_mG)
                delay(self.t_step_ms * ms)

        else:
            print("not yet implemented")

        ### finally, in case the worker refuses to die
        self.write_results()
