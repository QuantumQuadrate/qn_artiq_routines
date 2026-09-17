from collections import namedtuple

from artiq.experiment import *


Variable = namedtuple("Variable", ["name", "value", "value_type", "kwargs", "group"])


MASTER_SATELLITE_VARIABLES = (
    Variable(
        "n_measurements",
        100,
        NumberValue,
        {"type": "int", "ndecimals": 0, "step": 1, "scale": 1},
        "general",
    ),
    Variable(
        "t_delay_in_bob_mu",
        189,
        NumberValue,
        {"type": "int", "ndecimals": 0, "step": 1, "scale": 1},
        "master-satellite timing",
    ),
    Variable(
        "parallel_AOM_feedback",
        True,
        BooleanValue,
        {},
        "master-satellite feedback",
    ),
    # Two-atom thresholds are GLOBAL on purpose: they describe the joint
    # two-node readout, not one node's detector, so both nodes must agree on
    # them. They were per-node until 2026-09-17 and carried identical values
    # on both, which is the duplication this removes. Being global also makes
    # the bare name exist in two_nodes mode -- the per-node declarations only
    # ever produced two_atom_threshold_NodeX, while the two-node experiment
    # functions read self.two_atom_threshold.
    Variable(
        "two_atom_threshold",
        39000.0,
        NumberValue,
        {"type": "float"},
        "Thresholds and cut-offs",
    ),
    Variable(
        "two_atom_threshold_for_loading",
        89000.0,
        NumberValue,
        {"type": "float"},
        "Thresholds and cut-offs",
    ),
)


class ExperimentVariablesMasterSatelliteGlobal(EnvExperiment):
    """ExperimentVariables_master_satellite_global

    Initialize or intentionally update shared master-satellite variables.
    """

    def build(self):
        self.vars_list = list(MASTER_SATELLITE_VARIABLES)
        for variable in self.vars_list:
            try:
                value = self.get_dataset(variable.name)
            except KeyError:
                value = variable.value
            self.setattr_argument(
                variable.name,
                variable.value_type(value, **variable.kwargs),
                variable.group,
            )

    def run(self):
        for variable in self.vars_list:
            self.set_dataset(
                variable.name,
                getattr(self, variable.name),
                broadcast=True,
                persist=True,
            )
