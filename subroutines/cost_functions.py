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
    :param self: experiment instance
    :param photocounts: sequence containing photon count values for the measurement interval
    :return: -1 * atoms_loaded, the negated number of atoms loaded
    """

    atoms_loaded = 0
    q_last = (self.photocounts[0] > self.atom_counts_threshold)
    for x in self.photocounts[1:]:
        q = x > self.atom_counts_threshold
        if q != q_last and q_last:
            atoms_loaded += 1
        q_last = q
    atoms_loaded += q_last
    return -1 * atoms_loaded

def atom_loading_rate_pulsed_MOT_cost(self) -> TFloat:
    """
    the cost function for optimizing atom loading in a pulsed-MOT experiment

    by pulsed-MOT experiment, I mean the standard single atom experiment format which loads
    the dipole trap from the MOT, then turns off the MOT and checks for an atom with
    one readout, optionally doing more things following this before a second readout. I.e.,
    this is not for optimizing single atom loading in a steady state MOT

    :param self: experiment instance
    :return: -1*atom_retention, the negated number of atoms detected in the readout
    """

    shot1 = self.counts1_list
    atoms_loaded = [x > self.atom_counts_threshold for x in shot1]
    n_atoms_loaded = sum(atoms_loaded)
    loading_fraction = n_atoms_loaded / self.n_measurement
    return -100 * loading_fraction


def atom_retention_cost(self) -> TFloat:
    """
    the cost function for optimizing the fraction of atoms retained in a two-shot experiment.

    :param self: experiment instance
    :return: -1*retention_fraction, the negated number of atoms detected in the readout
    """

    shot1 = self.counts_list
    shot2 = self.counts2_list
    atoms_loaded = [x > self.atom_counts_threshold for x in shot1]
    n_atoms_loaded = sum(atoms_loaded)
    atoms_retained = [x > self.atom_counts2_threshold and y for x, y in zip(shot2, atoms_loaded)]
    retention_fraction = 0 if not n_atoms_loaded > 0 else sum(atoms_retained) / n_atoms_loaded

    return -100 * retention_fraction