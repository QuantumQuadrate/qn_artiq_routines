"""Microwave transition scans and optimization for the master-satellite stack.

Ported from the standalone MicrowaveScanOptimizer, which stays untouched and
keeps serving the standalone systems (including the standalone HealthCheck
GVS files, whose expids point at the standalone class). scan_dict,
scan_options, and fit_model_dict are shared with the standalone module so the
scan definitions cannot drift.

which_node selects Node1 or Node2 (two_nodes is deliberately not supported
and never will be). The reused single-node experiment functions run through
the ordinary legacy projection; every persisted calibration - the resonance
centers, pi-pulse times, and the health_check_uw_* fidelities - is stored
node-suffixed through the dataset-redirect layer, so Node1 and Node2 results
stay separate. fit_parameter_* and parent_rid remain unsuffixed run state.

Resubmission contract: this experiment resubmits itself (health-check
failures) with its own file/class name and carries which_node forward.
Future master-satellite HealthCheck files must do the same; the standalone
HealthCheck files keep targeting the standalone optimizer.
"""

import ast
from pathlib import Path

from artiq.experiment import *
from artiq.language import us, ns, MHz
import logging

import numpy as np
from numpy import array  # necessary for some override_ExperimentVariable entries

import copy

# this is where your experiment function should live
from subroutines.experiment_functions import *
import subroutines.experiment_functions as exp_functions

#### fitting functions (same models as the standalone optimizer)
from fitting.rabi_flop import rabi_flop
from fitting.rabi_flop_reversed import rabi_flop_reversed
from fitting.resonance_dip import resonance_dip
from fitting.resonance_peak import resonance_peak

from utilities.BaseExperiment_master_satellite import (
    BaseExperimentMasterSatellite,
    merge_per_node_override_dictionaries,
    _DatasetRedirectMixin,
)


fit_model_dict = {
    "rabi_flop": rabi_flop,
    "rabi_flop_reversed": rabi_flop_reversed,
    "resonance_dip": resonance_dip,
    "resonance_peak": resonance_peak,
}


def _load_standalone_scan_definitions():
    """Read scan_options/scan_dict from the standalone optimizer's source.

    The standalone MicrowaveScanOptimizer.py stays the single source of the
    scan definitions so the two versions cannot drift, but importing that
    module would drag the whole standalone Base/write_h5 stack into every
    repository-scan worker. Both assignments are pure literals, so they are
    extracted from the source instead.
    """
    source = (
        Path(__file__).parent / "standalone" / "MicrowaveScanOptimizer.py"
    ).read_text(encoding="utf-8")
    definitions = {}
    for node in ast.parse(source).body:
        if (
            isinstance(node, ast.Assign)
            and isinstance(node.targets[0], ast.Name)
            and node.targets[0].id in ("scan_options", "scan_dict")
        ):
            definitions[node.targets[0].id] = ast.literal_eval(node.value)
    return definitions["scan_options"], definitions["scan_dict"]


scan_options, scan_dict = _load_standalone_scan_definitions()


class MicrowaveScanOptimizer_master_satellite(
    _DatasetRedirectMixin, EnvExperiment
):
    """MicrowaveScanOptimizer_master_satellite

    Scan and optimize microwave transitions on selected Node1 or Node2.
    """

    VALID_NODES = ("Node1", "Node2")
    # Tests inject a fake stabilizer class here; None selects the real
    # AOMPowerStabilizer from subroutines/aom_feedback.py.
    _stabilizer_factory = None

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        # Repository examination supplies None for arguments and may run with
        # an empty dataset namespace; bind the suffixed device superset now
        # and defer mode validation and dataset loading to prepare().
        self.base = BaseExperimentMasterSatellite(experiment=self)
        self.base.build()

        self.setattr_argument(
            "which_node", EnumerationValue(self.VALID_NODES), "Node selection"
        )

        ## For Bookeeping purpose
        self.setattr_argument("parent_rid", NumberValue(0, ndecimals=0, step=1), "** For bookeeping purpose - always set this to 0 **")

        self.setattr_argument('run_health_check_and_optimize', BooleanValue(default=True), "Health Check")
        self.setattr_argument("target_fidelity", NumberValue(0.80, ndecimals=2, step=1), "Health Check")
        self.setattr_argument(
            'show_shared_applets',
            BooleanValue(default=True),
            "Applets",
        )

        self.setattr_argument("n_measurements", NumberValue(100, ndecimals=0, step=1), "General Scan Setting")
        # One override dictionary per node, and deliberately no shared one:
        # both may stay populated at once, so switching which_node is the only
        # edit needed to move a scan between nodes. Globals such as
        # n_measurements go in the running node's dictionary.
        for node in self.VALID_NODES:
            self.setattr_argument(
                f'override_ExperimentVariables_{node}',
                StringValue("{}"),
                "General Scan Setting",
                tooltip=f"Overrides applied only when {node} runs. Otherwise "
                        f"ignored entirely -- not even evaluated, so it may "
                        f"name {node} variables that do not exist on the "
                        f"other node. Globals belong here too.",
            )
        self.setattr_argument('enable_fitting', BooleanValue(default=True), "General Scan Setting")
        self.setattr_argument('enable_geometric_frequency_scan', BooleanValue(default=True), "General Scan Setting")

        for scan_name in scan_options:
            self.setattr_argument(
                scan_name,
                BooleanValue(default=False),
                "Microwave Scans - Select one of the following",
            )

        #### Frequency scan variables
        self.setattr_argument("freq_scan_range_left_kHz", NumberValue(125.0, ndecimals=1, step=1), "[Default] Frequency Scan variables - centered at resonance")
        self.setattr_argument("freq_scan_range_right_kHz", NumberValue(125.0, ndecimals=1, step=1), "[Default] Frequency Scan variables - centered at resonance")
        self.setattr_argument("freq_scan_step_size_kHz", NumberValue(20.0, ndecimals=1, step=1), "[Default] Frequency Scan variables - centered at resonance")

        #### Frequency scan variables - for faster scan!
        self.setattr_argument("shrink_factor", NumberValue(2.5, ndecimals=1, step=1), "[Geometric] Frequency Scan variables - centered at resonance")
        self.setattr_argument("freq_scan_half_range_kHz", NumberValue(125.0, ndecimals=1, step=1), "[Geometric] Frequency Scan variables - centered at resonance")
        self.setattr_argument("freq_scan_min_step_size_kHz", NumberValue(10.0, ndecimals=1, step=1), "[Geometric] Frequency Scan variables - centered at resonance")

        #### Time scan variables
        self.setattr_argument("time_scan_sequence", StringValue('np.arange(0,10,1)*us'), "Time Scan variables")

    def prepare(self):
        """
        performs initial calculations and sets parameter values before
        running the experiment.
        """
        # n_measurements is a run-local GUI value sharing its name with the
        # persistent global dataset; capture it before configure_execution
        # loads that dataset over the same attribute.
        self._execution_n_measurements = self.n_measurements

        selected_node = str(self.which_node)
        if selected_node not in self.VALID_NODES:
            raise ValueError(
                f"Unsupported which_node {self.which_node!r}; expected "
                "'Node1' or 'Node2'."
            )
        self._selected_node = selected_node

        self.base.configure_execution("single_node", selected_node)
        self.n_measurements = self._execution_n_measurements

        # The reused experiment functions branch on the alice/bob
        # presentation.
        self.which_node = self.base.NODE_LEGACY_NAMES[selected_node]
        self.base.prepare()
        self.base.prepare_laser_stabilizer(
            stabilizer_factory=self._stabilizer_factory
        )

        ### If this is the very first run (user-submitted), parent_rid will be 0
        if int(self.parent_rid) == 0:
            self.parent_rid = int(self.scheduler.rid)
            print("[Parent initialized] This run is the parent. Setting parent_rid =", self.parent_rid)
        else:
            print("[Child scan] Inheriting parent_rid =", self.parent_rid)

        # parent_rid is run bookkeeping and deliberately stays unsuffixed.
        self.set_dataset("parent_rid", self.parent_rid, broadcast=True, persist=True)

        ### goes through the booleans and sets the scan type
        self.scan_type = self.get_scan_type()

        # User overrides carry literal values; the scan definition's
        # override_items name a source attribute (e.g. t_microwave_00_pulse)
        # that is re-read at every application, matching the standalone
        # health-check re-merge behavior. Overrides mutate the suffixed
        # authoritative attributes; the legacy projection then mirrors them
        # for the reused functions.
        #
        # Only the running node's dictionary is read, so the other node's may
        # stay populated. The merge is shared with the GVS mixin through the
        # base module; globals() keeps override text evaluating in THIS
        # module's namespace, where entries like {'f_FORT': 240*MHz} find
        # their units.
        self._user_authoritative_overrides = (
            merge_per_node_override_dictionaries(
                experiment=self,
                active_nodes=(self._selected_node,),
                resolve_target=self.base.resolve_experiment_variable_target,
                eval_globals=globals(),
            )
        )
        self._scan_override_items = dict(
            scan_dict[self.scan_type]["override_items"]
        )
        self._apply_run_wide_overrides()

        ### setting the scan variable: legacy name for labels/metadata, the
        ### suffixed authoritative attribute for the actual scan writes.
        self.scan_variable1_name = scan_dict[self.scan_type]["scan_variable1_name"]
        self.scan_variable1 = self.base.resolve_experiment_variable_target(
            str(self.scan_variable1_name)
        )

        ### setting the experiment function (same for both nodes)
        self.experiment_name = scan_dict[self.scan_type]["experiment_name"]
        self._selected_experiment_function = getattr(
            exp_functions, self.experiment_name
        )
        self.experiment_function = (
            lambda: self._selected_experiment_function(self)
        )

        ### setting up the fit model
        fit_model = scan_dict[self.scan_type]["fit_model"]
        self.which_fit_model = fit_model_dict[fit_model]
        self.initialise = scan_dict[self.scan_type]["initialise"]

        ### setting up the scan sequence
        if self.scan_type.startswith('Freq'):
            center = getattr(self, scan_dict[self.scan_type]["center"])
            self.center = center
            print("center: ", center)
            if self.enable_geometric_frequency_scan:
                self.scan_sequence1 = self.make_scan_list(center=center,
                                                          half_range_kHz=self.freq_scan_half_range_kHz,
                                                          min_step_kHz=self.freq_scan_min_step_size_kHz,
                                                          mode="pair")
            else:
                self.scan_sequence1 = np.arange(center - self.freq_scan_range_left_kHz * kHz,
                                                center + self.freq_scan_range_right_kHz * kHz,
                                                self.freq_scan_step_size_kHz * kHz)

        elif self.scan_type.startswith(('Time', 'Ramsey')):
            self.scan_sequence1 = eval(self.time_scan_sequence)
        else:
            raise KeyError("There is no such scan type: ", self.scan_type)

        ### Note: some lines exist so the results stay compatible with the
        ### GVS analysis code.
        self.n_iterations1 = len(self.scan_sequence1)
        self.scan_variable2_name = ''
        self.scan_variable2 = None
        self.scan_sequence2 = np.zeros(1)
        self.n_iterations2 = 0

        scan_vars = [self.scan_variable1_name, self.scan_variable2_name]
        scan_vars = [x for x in scan_vars if x != '']
        self.scan_var_labels = ','.join(scan_vars)
        self.scan_var_filesuffix = '_and_'.join(scan_vars)

        self.set_dataset('scan_variable1_name', self.scan_variable1_name, broadcast=True)
        self.set_dataset('scan_sequence1', self.scan_sequence1, broadcast=True)
        self.set_dataset('scan_variable2_name', self.scan_variable2_name, broadcast=True)
        self.set_dataset('scan_sequence2', self.scan_sequence2, broadcast=True)
        self.set_dataset('experiment_function', self.experiment_name, broadcast=True)

        self.measurement = 0
        self.SPCM0_RO1 = 0
        self.SPCM0_RO2 = 0
        self.SPCM1_RO1 = 0
        self.SPCM1_RO2 = 0

        self.SPCM0_OtherNode_RO1 = 0
        self.SPCM0_OtherNode_RO2 = 0
        self.SPCM1_OtherNode_RO1 = 0
        self.SPCM1_OtherNode_RO2 = 0

        # If earlier queued experiments may still update datasets, reload the
        # authoritative values at run start (the master-satellite Base
        # prepare is one-shot, unlike the standalone rebuild pattern).
        status_dict = self.scheduler.get_status()
        my_rid = self.scheduler.rid
        earlier_experiments = len([rid for rid, _ in status_dict.items() if rid < my_rid])
        logging.info(
            "RID %s has %s earlier queued experiment(s)",
            my_rid, earlier_experiments,
        )
        self.needs_fresh_build = earlier_experiments > 0

    def _refresh_runtime_variable_state(self):
        self.base.refresh_compatibility_variables()
        self.base.refresh_variable_dependent_state()

    def _apply_run_wide_overrides(self):
        overrides = dict(self._user_authoritative_overrides)
        for legacy_name, source_name in self._scan_override_items.items():
            overrides[
                self.base.resolve_experiment_variable_target(str(legacy_name))
            ] = getattr(self, source_name)

        self.authoritative_overrides = overrides
        for target, value in overrides.items():
            setattr(self, target, value)
        self._refresh_runtime_variable_state()

    def _apply_scan_point(self, variable1_value):
        setattr(self, self.scan_variable1, variable1_value)
        self._refresh_runtime_variable_state()

    @kernel
    def initialize_hardware(self):
        self.core.break_realtime()
        self.base.initialize_hardware()

    @kernel
    def force_other_node_off(self):
        self.base.force_other_node_off()

    def _create_shared_applets(self):
        """Keep the node-independent applets up, if enabled.

        Microwaves Health Check and feedback RF, one of each per node. This
        experiment deliberately does not take over the per-node applet set or
        retire the idle node's group, which is what create_applets_for would
        do -- it is these plots we want on screen, not the whole dashboard.

        Never fatal: an applet is a convenience, and the CCB is a
        dashboard-side service that artiq_run and the offline compile check
        do not provide meaningfully. In particular, nothing is created when
        no dashboard is connected.
        """
        if not self.show_shared_applets:
            return
        try:
            from applets_master_satellite import (
                SCAN_APPLET_SPECS,
                SHARED_APPLET_SPECS,
                create_shared_applets_for,
            )

            # SCAN_APPLET_SPECS adds "Time without atom (s)". It rides the
            # shared path rather than the per-node one so this experiment and
            # GVS create one applet between them, not one each under different
            # groups -- applet identity is (name, group).
            created = create_shared_applets_for(
                self, self.base,
                shared_specs=SHARED_APPLET_SPECS + SCAN_APPLET_SPECS,
            )
        except Exception as error:  # noqa: BLE001 - convenience only
            logging.warning("could not create shared applets: %s", error)
        else:
            logging.info("requested %d shared applets", len(created))

    def initialize_datasets(self):
        self.base.initialize_result_state()

        self.set_dataset(self.scan_var_dataset, self.scan_var_labels, broadcast=True)
        self.set_dataset(self.scan_sequence1_dataset, self.scan_sequence1, broadcast=True)
        self.set_dataset(self.scan_sequence2_dataset, self.scan_sequence2, broadcast=True)

        # Archive the applied overrides under their authoritative names.
        for var, val in self.authoritative_overrides.items():
            self.set_dataset(var, val)

        # The stabilizer dB histories were already seeded node-suffixed by
        # Base.prepare_laser_stabilizer().

    def reset_datasets(self):
        """
        set datasets that are redefined each iteration.

        typically these datasets are used for plotting which would be
        meaningless if we continued to append to the data.
        """
        self.set_dataset("test_dataset", [0], broadcast=True)
        self.set_dataset('SPCM0_RO1_current_iteration', [0], broadcast=True)
        self.set_dataset('SPCM0_RO2_current_iteration', [0], broadcast=True)
        self.set_dataset('SPCM1_RO1_current_iteration', [0], broadcast=True)
        self.set_dataset('SPCM1_RO2_current_iteration', [0], broadcast=True)

        self.set_dataset('SPCM0_OtherNode_RO1_current_iteration', [0], broadcast=True)
        self.set_dataset('SPCM0_OtherNode_RO2_current_iteration', [0], broadcast=True)
        self.set_dataset('SPCM1_OtherNode_RO1_current_iteration', [0], broadcast=True)
        self.set_dataset('SPCM1_OtherNode_RO2_current_iteration', [0], broadcast=True)

        self.set_dataset('AllSPCMs_RO1_current_iteration', [0], broadcast=True)
        self.set_dataset('AllSPCMs_RO2_current_iteration', [0], broadcast=True)

        self.set_dataset('AllSPCMs_alternating_RO_alice_current_iteration', [0], broadcast=True)
        self.set_dataset('AllSPCMs_alternating_RO_bob_current_iteration', [0], broadcast=True)

        # these are set here because hardware initialization resets these
        self.set_dataset(self.scan_var_dataset, self.scan_var_labels, broadcast=True)
        self.set_dataset(self.scan_sequence1_dataset, self.scan_sequence1, broadcast=True)
        self.set_dataset(self.scan_sequence2_dataset, self.scan_sequence2, broadcast=True)

    def run(self):
        """
        Step through the scan sequence and run the experiment function,
        re-initializing hardware each iteration.
        """
        if self.needs_fresh_build:
            # Reload the authoritative datasets that earlier queued
            # experiments may have updated; keep the run-local GUI
            # n_measurements and re-apply the run-wide overrides.
            self.base.reload_experiment_variables()
            self.n_measurements = self._execution_n_measurements
            self._apply_run_wide_overrides()

        self.initialize_datasets()

        # Re-assert the two Microwaves Health Check applets (one per node).
        # They sit in their own group and are never retired, so this only has
        # to ask for them; they persist between runs and across node changes.
        self._create_shared_applets()

        # Once per run, ahead of both the health-check and the scan branches:
        # single_node mode leaves the other node out of the initialization
        # lifecycle, so without this its beams and coils keep whatever state
        # the previous run left them in. The per-scan-point
        # initialize_hardware() below cannot undo it -- core.reset() only ever
        # drives TTLs low and cannot change a held Zotino output.
        self.force_other_node_off()

        iteration = 0
        self.set_dataset("iteration", iteration, broadcast=True)

        if self.scan_type.startswith("Freq") and self.run_health_check_and_optimize:
            if self.health_check_general() == False:
                # write and overwrite the health check results
                self.write_results({'name': self.base.result_name_tag() + "_parent_rid_" + f"{self.parent_rid}" + "_" + self.experiment_name[
                                                                                         :-11] + "_scan_over_" + self.scan_var_filesuffix})

                print("Initial Health Check - failed with fidelity: ",
                      self.get_dataset(scan_dict[self.scan_type]["health_check_dataset_name"]))
                print("Scheduling Experiment for Optimization...")
                self.submit_opt_exp_general()

                self.scheduler.request_termination(self.scheduler.rid)

            self.core.comm.close()  # placing the hardware in a safe state and disconnecting it from the core device
            self.scheduler.pause()

        else:
            ##### Scan sequence - same as GVS (compatible with GVS analysis)
            if self.scan_type.startswith("Freq") and self.enable_geometric_frequency_scan:
                scan_with_fixed_sequence = False    ### adaptive scan enabled if False
            else:
                scan_with_fixed_sequence = True

            if scan_with_fixed_sequence:
                print("scanning with fixed sequence")
                for variable1_value in self.scan_sequence1:
                    self._apply_scan_point(variable1_value)
                    logging.info(f"current iteration: {self.scan_variable1_name} = {variable1_value}")

                    self.set_dataset("iteration", iteration, broadcast=True)

                    self.initialize_hardware()
                    self.reset_datasets()

                    # the measurement loop.
                    self.experiment_function()

                    # write and overwrite the file here so we can quit the experiment early without losing data
                    self.write_results({'name': self.base.result_name_tag() + "_parent_rid_" + f"{self.parent_rid}" + "_" + self.experiment_name[
                                                                                             :-11] + "_scan_over_" + self.scan_var_filesuffix})

                    iteration += 1

            else:
                print("scanning with adpative sequence")
                pending_sequence1 = list(self.scan_sequence1)
                did_refine = False

                scanned_points = []  # store every scanned variable1 value

                while pending_sequence1:
                    variable1_value = pending_sequence1.pop(0)

                    # ----- inner scanning loop -----
                    self._apply_scan_point(variable1_value)

                    logging.info(f"current iteration: {self.scan_variable1_name} = {variable1_value}")

                    self.set_dataset("iteration", iteration, broadcast=True)

                    self.initialize_hardware()
                    self.reset_datasets()

                    self.experiment_function()
                    self.write_results({'name': self.base.result_name_tag() + "_parent_rid_" + f"{self.parent_rid}" + "_" + self.experiment_name[
                                                                                             :-11] + "_scan_over_" + self.scan_var_filesuffix})

                    iteration += 1

                    # store the scanned variable1_value
                    scanned_points.append(variable1_value)

                    self.scan_sequence1 = np.array(scanned_points)
                    self.set_dataset("scan_sequence1", (scanned_points), broadcast=True)

                    if (not did_refine) and (len(scanned_points) == 4):
                        print("in decision loop")

                        AllSPCMs_RO1 = self.get_dataset("AllSPCMs_RO1")
                        AllSPCMs_RO2 = self.get_dataset("AllSPCMs_RO2")

                        retention_array, loading_rate_array, n_atoms_loaded_array = self.get_loading_and_retention(
                            AllSPCMs_RO1, AllSPCMs_RO2,
                            self.n_measurements,
                            int((len(AllSPCMs_RO1) - 1) / self.n_measurements),
                            self.single_atom_threshold * self.t_SPCM_first_shot
                        )

                        fit_model = scan_dict[self.scan_type]["fit_model"]
                        ### setting target retention to calculate fidelity
                        if fit_model in ["resonance_dip", "rabi_flop"]:
                            pass
                        elif fit_model in ["resonance_peak", "rabi_flop_reversed"]:
                            retention_array = 1.0 - retention_array

                        # First four points of the PAIRED scan sequence
                        # (make_scan_list, mode="pair"):
                        #   index 0 -> -x, 1 -> +x, 2 -> -r*x, 3 -> +r*x
                        # where r = 1 - 1/shrink_factor. These are NOT at
                        # -x/2 and +x/2 unless shrink_factor is exactly 2; at
                        # shrink_factor 2.5 they sit at +-0.6x, so any
                        # threshold tuned on the "half" reading is off by 20%.
                        threshold_low = 0.7  # tune this as needed
                        new_center = None
                        no_clear_dip = False
                        # NaN so the all(> 0.9) broaden test below is False if
                        # we never got four points.
                        R0 = R1 = R2 = R3 = float("nan")

                        if len(retention_array) < 4 or len(self.scan_sequence1) < 4:
                            # A run can end with fewer data windows than scan
                            # points, so four are not guaranteed.
                            print(
                                "Decision loop: only "
                                f"{min(len(retention_array), len(self.scan_sequence1))}"
                                " point(s) available; leaving the center alone"
                            )
                        else:
                            f0, f1, f2, f3 = self.scan_sequence1[0:4]
                            R0, R1, R2, R3 = map(float, retention_array[0:4])

                            # WING ASYMMETRY, not slope, locates the
                            # resonance. The old test (slope_left < 0 and
                            # slope_right > 0) only says "some dip lies
                            # between the outer points", which holds for any
                            # interior dip however lopsided -- that is why a
                            # clearly right-shifted scan was reported as well
                            # centered. Positive asymmetry means the right
                            # probe is deeper, i.e. the resonance is RIGHT of
                            # the current center.
                            asym_inner = R2 - R3
                            asym_outer = R0 - R1

                            # Retention is a fraction over ~n_measurements
                            # loaded atoms, so its binomial sigma is roughly
                            # sqrt(p(1-p)/n) ~ 0.05 at n = 100, and ~0.07 on
                            # a difference of two. This is about a 2-sigma
                            # cut: raise it if the decision looks jumpy,
                            # lower it to react to smaller shifts.
                            asym_tol = 0.15

                            # Decide on whichever PAIR actually shows the
                            # asymmetry. The inner pair alone is not enough: a
                            # resonance at or beyond an outer probe leaves
                            # BOTH inner points on the flat shoulder, so
                            # asym_inner is ~0 while asym_outer is large. A
                            # dip sitting exactly on f0 (R = [0, .75, .75,
                            # .75]) reads as perfectly centered on the inner
                            # pair alone.
                            if abs(asym_outer) > abs(asym_inner):
                                asym = asym_outer
                            else:
                                asym = asym_inner

                            lowest_index = int(np.argmin(retention_array[:4]))
                            lowest_value = float(retention_array[lowest_index])

                            if abs(asym_outer) > asym_tol and asym_outer * asym_inner < 0:
                                # The inner and outer pairs disagree about
                                # which side is deeper. Usually means the dip
                                # is comparable to the noise, or the window
                                # straddles two features. Worth seeing.
                                print(
                                    "Decision loop: inner and outer wings "
                                    f"disagree (inner {asym_inner:+.3f}, "
                                    f"outer {asym_outer:+.3f})"
                                )

                            # Skew is tested BEFORE centering. The centered
                            # test used to come first with no depth or
                            # symmetry guard, so it absorbed nearly every
                            # skewed scan and left cases 2.3/2.4 unreachable.
                            # Depth deliberately does not gate the skew cases:
                            # a resonance outside the probed window leaves all
                            # four points high, which is exactly when
                            # re-centering matters most.
                            if asym > asym_tol:
                                if lowest_index == 1:
                                    print("Case 2.2: right-skewed (out-of-range)")
                                    if lowest_value < 0.2:
                                        # freq_scan_min_step_size_kHz is a kHz
                                        # number; scan values are in Hz
                                        new_center = f1 + self.freq_scan_min_step_size_kHz * kHz
                                    elif lowest_value < 0.4:
                                        step = abs(f1 - f3)
                                        new_center = f1 + step
                                    else:
                                        step = 2 * abs(f1 - f3)
                                        new_center = f1 + step
                                else:
                                    print("Case 2.4: right-skewed (within-range)")
                                    if lowest_value < 0.2:
                                        new_center = f3 - self.freq_scan_min_step_size_kHz * kHz
                                    else:
                                        step = abs((f1 - f3) / 2)
                                        new_center = f3 - step

                            elif asym < -asym_tol:
                                if lowest_index == 0:
                                    print("Case 2.1: left-skewed (out-of-range)")
                                    if lowest_value < 0.2:
                                        new_center = f0 - self.freq_scan_min_step_size_kHz * kHz
                                    elif lowest_value < 0.4:
                                        step = abs(f2 - f0)
                                        new_center = f0 - step
                                    else:
                                        step = 2 * abs(f2 - f0)
                                        new_center = f0 - step
                                else:
                                    print("Case 2.3: left-skewed (within-range)")
                                    if lowest_value < 0.2:
                                        new_center = f2 + self.freq_scan_min_step_size_kHz * kHz
                                    else:
                                        step = abs((f2 - f0) / 2)
                                        new_center = f2 + step

                            elif lowest_value > threshold_low:
                                # Symmetric AND shallow: no usable dip in the
                                # probed window at all. Hands over to the
                                # broaden-the-search path below rather than
                                # silently accepting the center.
                                no_clear_dip = True
                                print(
                                    "Decision loop: no clear dip in the probed "
                                    f"window (best retention {lowest_value:.3f} "
                                    f"> {threshold_low})"
                                )

                            else:
                                print(
                                    "Case 2.0: Well centered (asymmetry "
                                    f"{asym_inner:+.3f} within {asym_tol})"
                                )

                        # If we decided on a refined center, build a new scan
                        if new_center is not None:
                            print(f"Changing center to {new_center}")
                            new_sequence1 = self.make_scan_list(
                                center=new_center,
                                half_range_kHz=self.freq_scan_half_range_kHz,
                                min_step_kHz=self.freq_scan_step_size_kHz,
                                mode="pair",
                            )
                            pending_sequence1 = list(new_sequence1)
                            did_refine = True
                        else:
                            if no_clear_dip or all(r > 0.9 for r in (R0, R1, R2, R3)):
                                print("Resonance not in this range; searching broader range with same center")
                                self.shrink_factor = self.shrink_factor * 2
                                print("Manual scan range: - 2* self.freq_scan_range_left_kHz ~ +2* self.freq_scan_range_left_kHz, 20kHz step")

                                new_sequence1 = np.arange(self.center - 2 * self.freq_scan_range_left_kHz * kHz,
                                                          self.center + 2 * self.freq_scan_range_right_kHz * kHz,
                                                          20 * kHz)

                                pending_sequence1 = list(new_sequence1)
                                did_refine = True

            #### for fitting
            if self.enable_fitting and not self.scan_type.startswith("Ramsey"):
                AllSPCMs_RO1 = self.get_dataset("AllSPCMs_RO1")
                AllSPCMs_RO2 = self.get_dataset("AllSPCMs_RO2")
                retention_array = self.get_retention(AllSPCMs_RO1, AllSPCMs_RO2, self.n_measurements, len(self.scan_sequence1),
                                   self.single_atom_threshold * self.t_SPCM_first_shot)

                t = self.scan_sequence1
                y = retention_array

                p, p_err = self.which_fit_model.fit(t, y, evaluate_function=False, initialise=self.initialise)
                print(p)

                ### saving fitted parameters (unsuffixed run state, per plan)
                for var, val in p.items():
                    ds_name = f"fit_parameter_{var}"
                    self.set_dataset(ds_name, round(float(val), 7), broadcast=True, persist=True)

                for var, val in p_err.items():
                    ds_name = f"fit_parameter_{var}_err"
                    self.set_dataset(ds_name, round(float(val), 7), broadcast=True, persist=True)

                ### depending on scan_type, the fitted parameter is different
                if self.scan_type.startswith("Freq"):
                    if self.health_check_general(fit_check=True):
                        print("optimization success - dataset updating to fitted value")

                        # Stored node-suffixed through the dataset redirect.
                        self.set_dataset(scan_dict[self.scan_type]["center"], round(float(p["f0"]), -3), broadcast=True, persist=True)
                        print(scan_dict[self.scan_type]["center"], " updated to ", getattr(self, scan_dict[self.scan_type]["center"]))

                    else:
                        print("optimization failed - dataset not updated")

                elif self.scan_type.startswith("Time"):
                    if self.health_check_general(fit_check=True):
                        print("optimization success - dataset updating to fitted value")
                        self.set_dataset(scan_dict[self.scan_type]["pi_pulse"], round(float(p["t_pi"]), 7), broadcast=True,
                                         persist=True)
                        print(scan_dict[self.scan_type]["pi_pulse"], " updated to ",
                              getattr(self, scan_dict[self.scan_type]["pi_pulse"]))

                    else:
                        print("optimization failed - dataset not updated")

                # write and overwrite the file here so we can quit the experiment early without losing data
                self.write_results({'name': self.base.result_name_tag() + "_parent_rid_" + f"{self.parent_rid}" + "_" + self.experiment_name[:-11] + "_scan_over_" + self.scan_var_filesuffix})

                iteration += 1

        print("****************    General Variable Scan DONE   *****************")

    ################################### Health Check Functions ########################################

    def health_check_general(self, fit_check=False):
        """
        Perform a single-point health check on the microwave or time-domain
        scan; see the standalone MicrowaveScanOptimizer for the full
        description. Fidelities are stored in the node-suffixed
        health_check_uw_* datasets through the dataset redirect.
        """
        self.initialize_hardware()
        self.reset_datasets()

        if fit_check:
            if self.scan_type.startswith("Freq"):
                print("Original f0 value: ", getattr(self, scan_dict[self.scan_type]["center"]))
                fit_f0 = self.get_dataset("fit_parameter_f0")
                self._set_calibration_attribute(
                    scan_dict[self.scan_type]["center"], round(float(fit_f0), -3)
                )
                print("Fit check with fitted f0: ", getattr(self, scan_dict[self.scan_type]["center"]))
            elif self.scan_type.startswith("Time"):
                print("Original", scan_dict[self.scan_type]["pi_pulse"], "value: ", getattr(self, scan_dict[self.scan_type]["pi_pulse"]))
                fit_pi_pulse = self.get_dataset("fit_parameter_t_pi")
                self._set_calibration_attribute(
                    scan_dict[self.scan_type]["pi_pulse"], round(float(fit_pi_pulse), 7)
                )
                print("Fit check with fitted pi pulse: ", getattr(self, scan_dict[self.scan_type]["pi_pulse"]))
        else:
            self._apply_run_wide_overrides()

            if self.scan_type.startswith("Freq"):
                print("Health Check with original f0: ", getattr(self, scan_dict[self.scan_type]["center"]))
            elif self.scan_type.startswith("Time"):
                print("Health Check with original pi pulse: ", getattr(self, scan_dict[self.scan_type]["pi_pulse"]))

        if self.scan_type.startswith("Freq"):
            check_value = getattr(self, scan_dict[self.scan_type]["center"])
        elif self.scan_type.startswith("Time"):
            check_value = getattr(self, scan_dict[self.scan_type]["pi_pulse"])
        self._apply_scan_point(check_value)

        ### setting experiment function (same for both nodes)
        self.experiment_name = scan_dict[self.scan_type]["experiment_name"]

        ### setting target retention to calculate fidelity
        if scan_dict[self.scan_type]["fit_model"] in ["resonance_dip", "rabi_flop"]:
            target_retention_0 = True
        elif scan_dict[self.scan_type]["fit_model"] in ["resonance_peak", "rabi_flop_reversed"]:
            target_retention_0 = False

        self.initialise = scan_dict[self.scan_type]["initialise"]

        self._selected_experiment_function = getattr(
            exp_functions, self.experiment_name
        )
        self.experiment_function = (
            lambda: self._selected_experiment_function(self)
        )
        self.experiment_function()

        AllSPCMs_RO1 = self.get_dataset("AllSPCMs_RO1")
        AllSPCMs_RO2 = self.get_dataset("AllSPCMs_RO2")

        retention_array, loading_rate_array, n_atoms_loaded_array = self.get_loading_and_retention(AllSPCMs_RO1,
                                                                                                   AllSPCMs_RO2,
                                                                                                   self.n_measurements,
                                                                                                   int((len(AllSPCMs_RO1) - 1) / self.n_measurements),
                                                                                                   self.single_atom_threshold * self.t_SPCM_first_shot)

        if target_retention_0:
            fidelity = 1.0 - retention_array[-1]
        else:
            fidelity = retention_array[-1]

        ##### update health_check_dataset only with FREQ scans;
        if self.scan_type.startswith("Freq"):
            self.set_dataset(scan_dict[self.scan_type]["health_check_dataset_name"], round(float(fidelity), 2), broadcast=True, persist=True)

        if fidelity >= self.target_fidelity:
            print("Health Check - passed with fidelity: ", fidelity)

            #### update health_check dataset with TIME scans only when passed
            if self.scan_type.startswith("Time"):
                self.set_dataset(scan_dict[self.scan_type]["health_check_dataset_name"], round(float(fidelity), 2),
                                 broadcast=True, persist=True)
            return True
        else:
            return False

    def _set_calibration_attribute(self, legacy_name, value):
        """Set a calibration through its authoritative suffixed attribute."""
        resolved = self.base.resolve_experiment_variable_target(
            str(legacy_name)
        )
        setattr(self, resolved, value)
        self._refresh_runtime_variable_state()

    def submit_opt_exp_general(self, override_arguments=None):
        """
        Schedules an optimization run of THIS master-satellite class when a
        health check fails, carrying which_node forward. Future
        master-satellite HealthCheck files must submit the same expid shape.
        """
        print("submitting another experiment")
        ## 99 seems to be the highest priority that can be set to.

        default_expid = {
            "log_level": 30,
            "file": "qn_artiq_routines\\MicrowaveScanOptimizer_master_satellite.py",
            "class_name": "MicrowaveScanOptimizer_master_satellite",
            "arguments": {
                "which_node": self._selected_node,
                "parent_rid": self.parent_rid,
                "run_health_check_and_optimize": False,
                "target_fidelity": 0.80,

                "n_measurements": self._execution_n_measurements,
                "override_ExperimentVariables_Node1": "{}",
                "override_ExperimentVariables_Node2": "{}",
                "enable_fitting": True,
                "enable_geometric_frequency_scan": True,

                # Scan type toggles
                "Frequency_00_Scan": False,
                "Frequency_01_Scan": False,
                "Frequency_11_Scan": False,
                "Frequency_m10_Scan": False,
                "Frequency_m11_Scan": False,
                "Time_00_Scan": False,
                "Time_01_Scan": False,
                "Time_11_Scan": False,
                "Time_m10_Scan": False,
                "Time_m11_Scan": False,
                "Ramsey_00_Scan": False,
                "Ramsey_01_Scan": False,
                "Ramsey_11_Scan": False,

                # [Full Scan] Frequency scan parameters (kHz)
                "freq_scan_range_left_kHz": 100.0,
                "freq_scan_range_right_kHz": 101.0,
                "freq_scan_step_size_kHz": 20.0,

                # [Faster Scan] Frequency scan parameters (kHz)
                "shrink_factor": 2.5,
                "freq_scan_half_range_kHz": 100.0,
                "freq_scan_min_step_size_kHz": 10.0,

                # Time scan parameters
                "time_scan_sequence": 'np.arange(0,10,1)*us'
            },
            "repo_rev": "N/A",
        }

        new_expid = copy.deepcopy(default_expid)

        for key in scan_options:
            if getattr(self, key, False):     # if self.key does not exist, it returns False
                new_expid["arguments"][key] = True
                break  # stop once the first TRUE is found

        ### keeping BOTH per-node override dictionaries. A health-check
        ### resubmission must not drop the idle node's: it is carried by this
        ### expid alone, so losing it here would silently empty that field for
        ### every later run of the resubmitted chain.
        for node in self.VALID_NODES:
            argument_name = f"override_ExperimentVariables_{node}"
            new_expid["arguments"][argument_name] = getattr(self, argument_name)

        if override_arguments is not None:
            for key, value in override_arguments.items():
                if key in new_expid["arguments"]:
                    print(f"Overriding {key}: {new_expid['arguments'][key]} to {value}")
                    new_expid["arguments"][key] = value
                else:
                    raise KeyError(f"Invalid override key: '{key}'. ")

        self.scheduler.submit(pipeline_name="main", expid=new_expid, priority=99, due_date=None, flush=False)

    def make_scan_list(self, center, half_range_kHz, min_step_kHz=10.0, method="center_geometric", mode="pair"):
        """
        Create a geometric scan list around a center frequency; see the
        standalone MicrowaveScanOptimizer for the strategy description.
        """
        assert self.shrink_factor > 1, "shrink_factor must be > 1"
        assert mode in ("sequential", "pair"), f"Unknown mode: {mode}"

        x = float(half_range_kHz)
        r = (1 - 1 / self.shrink_factor)

        if method == "center_geometric":
            if mode == "sequential":
                offsets = []
                val = x
                prev = None
                while val >= min_step_kHz and (prev is None or abs(prev - val) >= min_step_kHz):
                    offsets.append(val)
                    prev = val
                    val *= r

                symmetric_offsets = sorted([-v for v in offsets] + [0.0] + offsets)
                values = [center + v * kHz for v in symmetric_offsets]

            elif mode == "pair":
                offsets = []
                val = x
                prev = None
                while val >= min_step_kHz and (prev is None or abs(prev - val) >= min_step_kHz):
                    offsets.append(-val)
                    offsets.append(+val)
                    prev = val
                    val *= r

                offsets.append(0.0)
                values = [center + v * kHz for v in offsets]

        elif method == "quarter_geometric":
            quarter = x / 2.0

            local_pos = []
            val_local = quarter
            while val_local >= min_step_kHz:
                local_pos.append(val_local)
                val_local *= r

            local_offsets = sorted([-v for v in local_pos] + [0.0] + local_pos)

            if mode == "sequential":
                raw_offsets = []
                for o in local_offsets:
                    raw_offsets.append(-quarter + o)
                for o in local_offsets:
                    raw_offsets.append(+quarter + o)
                raw_offsets.append(0.0)
                offsets = sorted({o for o in raw_offsets if -x <= o <= x})
                values = [center + v * kHz for v in offsets]

            elif mode == "pair":
                raw_offsets = [-quarter, +quarter]
                for o in local_offsets:
                    if o == 0.0:
                        continue
                    candidates = [
                        -quarter + o,
                        -quarter - o,
                        +quarter + o,
                        +quarter - o,
                    ]
                    for c in candidates:
                        if -x <= c <= x:
                            raw_offsets.append(c)

                if 0.0 not in raw_offsets:
                    raw_offsets.append(0.0)

                offsets = raw_offsets
                values = [center + v * kHz for v in offsets]

        return np.array(values)

    def get_scan_type(self):
        """
        Determine which scan type is enabled; exactly one flag must be True.
        """
        enabled = [name for name in scan_options if getattr(self, name, False)]

        if len(enabled) == 0:
            raise ValueError("No scan type selected. At least one scan flag must be True.")

        if len(enabled) > 1:
            raise ValueError(
                f"Multiple scan types selected ({enabled}). Only one can be True."
            )

        return enabled[0]

    def get_loading_and_retention(self, photocounts, photocounts2, measurements, iterations, cutoff):
        """
        Compute loading rate, retention, and number of atoms loaded per
        iteration.
        """
        retention_array = np.zeros(iterations)
        loading_rate_array = np.zeros(iterations)
        n_atoms_loaded_array = np.zeros(iterations)

        for i in range(iterations):
            # the SPCM datasets are seeded with a leading [0]; skip it so each
            # window holds exactly this iteration's measurements
            shot1 = photocounts[1 + i * measurements:1 + (i + 1) * measurements]
            shot2 = photocounts2[1 + i * measurements:1 + (i + 1) * measurements]

            atoms_loaded = [x > cutoff for x in shot1]
            n_atoms_loaded = sum(atoms_loaded)
            atoms_retained = [x > cutoff and y for x, y in zip(shot2, atoms_loaded)]
            retention_fraction = 0 if not n_atoms_loaded > 0 else sum(atoms_retained) / sum(atoms_loaded)
            loading_rate_array[i] = n_atoms_loaded / measurements
            n_atoms_loaded_array[i] = n_atoms_loaded
            retention_array[i] = retention_fraction
        return retention_array, loading_rate_array, n_atoms_loaded_array

    def get_retention(self, photocounts, photocounts2, measurements, iterations, cutoff):
        """
        Compute only the retention fraction per iteration.
        """
        retention_array = np.zeros(iterations)

        for i in range(iterations):
            # the SPCM datasets are seeded with a leading [0]; skip it so each
            # window holds exactly this iteration's measurements
            shot1 = photocounts[1 + i * measurements:1 + (i + 1) * measurements]
            shot2 = photocounts2[1 + i * measurements:1 + (i + 1) * measurements]

            atoms_loaded = [x > cutoff for x in shot1]
            n_atoms_loaded = sum(atoms_loaded)
            atoms_retained = [x > cutoff and y for x, y in zip(shot2, atoms_loaded)]
            retention_fraction = 0 if not n_atoms_loaded > 0 else sum(atoms_retained) / sum(atoms_loaded)
            retention_array[i] = retention_fraction

        return retention_array
