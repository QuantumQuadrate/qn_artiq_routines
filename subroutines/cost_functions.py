from artiq.experiment import *
from skimage.filters import threshold_otsu
import numpy as np

"""
Functions which can be used for optimization of various experiment variables

Note: for uniformity, I generally try to have the best cost be -100. This is
not required, but it makes it less likely that we will forget to correctly specify
the target cost to the optimizer

Adding new cost function:
1. The name of each function below MUST end in 'cost'
2. The cost function must take self as the first and only argument,
which is a reference to an experiment. This allows the cost function
to reference the self.SPCM0_RO1 we care about without needing to pass in self.SPCM0_RO1
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
    :param SPCM0_RO1: sequence containing photon count values for the measurement interval
    :return: -1 * atoms_loaded, the negated number of atoms loaded
    """

    atoms_loaded = 0
    q_last = (self.SPCM0_RO1[0] > self.single_atom_SPCM0_RO1_threshold)
    for x in self.SPCM0_RO1[1:]:
        q = x > self.single_atom_SPCM0_RO1_threshold
        if q != q_last and q_last:
            atoms_loaded += 1
        q_last = q
    atoms_loaded += q_last
    return -1 * atoms_loaded


def atom_loading_cost(self) -> TFloat:
    """
    the cost function for optimizing atom loading in a pulsed-MOT experiment

    by pulsed-MOT experiment, I mean the standard single atom experiment format which loads
    the dipole trap from the MOT, then turns off the MOT and checks for an atom with
    one readout, optionally doing more things following this before a second readout. I.e.,
    this is not for optimizing single atom loading in a steady state MOT

    :param self: experiment instance
    :return: -1*atom_retention, the negated number of atoms detected in the readout
    """

    shot1 = self.SPCM0_RO1_list
    atoms_loaded = [x > self.single_atom_SPCM0_RO1_threshold for x in shot1]
    n_atoms_loaded = sum(atoms_loaded)
    loading_fraction = n_atoms_loaded/len(shot1)

    cost = -100 * loading_fraction

    # todo delete
    print("inside cost function:", self.iteration, cost)
    return cost


def atom_loading_with_otsu_threshold_cost(self) -> TFloat:
    """
    the cost function for optimizing atom loading in a pulsed-MOT experiment

    by pulsed-MOT experiment, I mean the standard single atom experiment format which loads
    the dipole trap from the MOT, then turns off the MOT and checks for an atom with
    one readout, optionally doing more things following this before a second readout. I.e.,
    this is not for optimizing single atom loading in a steady state MOT

    :param self: experiment instance
    :return: -1*atom_retention, the negated number of atoms detected in the readout
    """

    shot1 = self.SPCM0_RO1_list
    atoms_loaded = [x > self.single_atom_SPCM0_RO1_threshold for x in shot1]
    n_atoms_loaded = sum(atoms_loaded)
    loading_fraction = n_atoms_loaded/len(shot1)

    # If there are enough atoms loaded according to the threshold, recompute the loading rate with an Otsu threshold.
    # this will typically give a more accurate cut-off in case the histogram cleanness or cut-off changes with the
    # parameters we are varying. The reason we can not start with Otsu thresholding right away is that it would
    # still return a cut-off even if we load no atoms, and the cut-off would just bisect the background mode. Put
    # another way, it can not tell whether the data is bimodal or not.
    if loading_fraction > 0.3:  # apparent very low rate loading might just be wrongly classified background
        threshold = threshold_otsu(np.array(self.SPCM0_RO1_list))
        atoms_loaded = [x > threshold for x in shot1]
        n_atoms_loaded = sum(atoms_loaded)
        loading_fraction = n_atoms_loaded / len(shot1)

    return -100 * loading_fraction


def atom_retention_and_loading_cost(self) -> TFloat:
    """
    the cost function for optimizing the fraction of atoms loaded and the fraction
    retained in a two-shot experiment.

    this cost ensures that retention is not optimized at the expense of loading,
    which is useful for tuning parameters that are relevant before or during the first readout

    :param self: experiment instance
    :return: -100*retention_fraction*loading_fraction/0.6
    """

    cost = 1
    shot1 = self.SPCM0_RO1_list
    shot2 = self.SPCM0_RO2_list
    atoms_loaded = [x > self.single_atom_SPCM0_RO1_threshold for x in shot1]
    n_atoms_loaded = sum(atoms_loaded)
    atoms_retained = [x > self.single_atom_SPCM0_RO2_threshold and y for x, y in zip(shot2, atoms_loaded)]
    retention_fraction = 0 if not n_atoms_loaded > 0 else sum(atoms_retained) / n_atoms_loaded
    loading_fraction = n_atoms_loaded/len(shot1)

    # recompute the retention and loading with an Otsu threshold.
    #todo: check if this correct

    # if loading_fraction > 0.3:  # apparent very low rate loading might just be wrongly classified background
    #     threshold = threshold_otsu(np.array(self.SPCM0_RO1_list))
    #     atoms_loaded = [x > threshold for x in shot1]
    #     n_atoms_loaded = sum(atoms_loaded)
    #     loading_fraction = n_atoms_loaded / len(shot1)
    #     retention_fraction = 0 if not n_atoms_loaded > 0 else sum(atoms_retained) / n_atoms_loaded
    #     cost *= threshold/500

    # 0.6 is probably the best loading rate we can hope for.
    # cost *= -50 * (retention_fraction + loading_fraction / 0.6)
    # now we are using atom_loading_2 which gives nearly 100% loading.
    cost *= -50 * (retention_fraction + loading_fraction)
    return cost


def atom_retention_cost(self) -> TFloat:
    """
    the cost function for optimizing the fraction of atoms retained in a two-shot experiment.

    :param self: experiment instance
    :return: -100*retention_fraction, the negated percentage of atoms detected in the readout
    """

    shot1 = self.BothSPCMs_RO1_list
    shot2 = self.BothSPCMs_RO2_list

    atoms_loaded = [x > self.single_atom_RO1_threshold for x in shot1]
    n_atoms_loaded = sum(atoms_loaded)
    atoms_retained = [x > self.single_atom_RO2_threshold and y for x, y in zip(shot2, atoms_loaded)]
    retention_fraction = 0 if not n_atoms_loaded > 0 else sum(atoms_retained) / n_atoms_loaded
    loading_fraction = n_atoms_loaded/len(shot1)

    # # recompute the retention and loading with an Otsu threshold.
    # # this will typically give a more accurate cut-off in case the histogram cleanness or cut-off changes with the
    # # parameters we are varying. The reason we can not start with Otsu thresholding right away is that it would
    # # still return a cut-off even if we load no atoms, and the cut-off would just bisect the background mode. Put
    # # another way, it can not tell whether the data is bimodal or not.
    # if loading_fraction > 0.3:  # apparent very low rate loading might just be wrongly classified background
    #     threshold = threshold_otsu(np.array(self.SPCM0_RO1_list))
    #     atoms_loaded = [x > threshold for x in shot1]
    #     n_atoms_loaded = sum(atoms_loaded)
    #     atoms_retained = [x > self.single_atom_SPCM0_RO2_threshold and y for x, y in zip(shot2, atoms_loaded)]
    #     retention_fraction = 0 if not n_atoms_loaded > 0 else sum(atoms_retained) / n_atoms_loaded

    return - 100 * retention_fraction
    # return -100 * (1-retention_fraction)
    ### use -100 * retention_fraction to maximize retention
    ### use -100 * (1 - retention_fraction) to minimize retention



def atom_blowaway_cost(self) -> TFloat:
    """
    the cost function for minimizing the fraction of atoms retained in a two-shot experiment.

    :param self: experiment instance
    :return: 100*(retention_fraction-1)
    """

    shot1 = self.SPCM0_RO1_list
    shot2 = self.SPCM0_RO2_list
    atoms_loaded = [x > self.single_atom_SPCM0_RO1_threshold for x in shot1]
    n_atoms_loaded = sum(atoms_loaded)
    atoms_retained = [x > self.single_atom_SPCM0_RO2_threshold and y for x, y in zip(shot2, atoms_loaded)]
    retention_fraction = 0 if not n_atoms_loaded > 0 else sum(atoms_retained) / n_atoms_loaded

    return 100*(retention_fraction - 1)

#
#
# def atom_lading_and_time_cost(self, data: TArray(TFloat,1)) -> TInt32:
#     atom_loading_time_list = self.Atom_loading_time
#
#     total_t = 0.0
#     for t in data:
#         total_t += -1.0 / t
#     average_t = total_t / len(data) ### Though I am naming these at _t, these are indeed 1/t to calculate the cost
#     return int(round(average_t))
#
#
# # use with this kind of format:  get_cost(self.atom_loading_time_list)