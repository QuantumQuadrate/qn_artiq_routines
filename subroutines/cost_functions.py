from artiq.experiment import *
import numpy as np
from abc import ABC

"""
Functions which can be used for optimization of various experiment variables

Adding new cost function:
1. The name of each function below MUST end in 'cost'
2. The cost function must take self as the first and only argument,
which is a reference to an experiment. This allows the cost function
to reference the self.counts we care about without needing to pass in self.counts
explicitly thus allowing us to keep the interface general.
3. The cost returned by the function must be a float.

See GeneralVariableOptimizer.py to use these functions.
"""

@kernel
def template_cost(self) -> TFloat:
    """
    """
    cost = 0.0

    return cost

def atoms_loaded_in_continuous_MOT_cost(self) -> TFloat:
    """
    the cost function for optimizing number of atoms in the dipole trap
    in a continuously loaded MOT
    :param self:
    :return:
    """

    q = 0 # whether the counts exceeded the atom threshold
    atoms_loaded = 0.0
    q_last = (self.counts[0] > self.atom_counts_threshold)
    for x in self.counts[1:]:
        q = x > self.atom_counts_threshold
        if q != q_last and q_last:
            atoms_loaded += 1
            # todo: replace with logging statement?
            # self.print_async("optimizer found the single atom signal!")
        q_last = q
    return -1.0 * atoms_loaded
