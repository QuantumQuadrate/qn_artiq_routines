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
)


class ExperimentVariablesMasterSatellite(EnvExperiment):
    """Initialize or intentionally update shared master-satellite variables."""

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
