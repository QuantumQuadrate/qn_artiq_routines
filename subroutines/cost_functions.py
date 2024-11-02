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


def atoms_loaded_in_continuous_MOT_cost(self) -> TInt32:
    """
    the cost function for optimizing number of atoms in the dipole trap
    in a continuously loaded MOT
    :param self: experiment instance
    :param photocounts: sequence containing photon count values for the measurement interval
    :return: -1 * atoms_loaded, the negated number of atoms loaded
    """
    data = self.counts_list
    atoms_loaded = 0
    q_last = (data[0] > self.atom_counts_threshold)
    for x in data[1:]:
        q = x > self.atom_counts_threshold
        if q != q_last and q_last:
            atoms_loaded += 1
        q_last = q
    atoms_loaded += q_last

    # recompute if a certain fraction of counts is above threshold.
    # fraction is somewhat arbitrary should be high enough to mitigate false positives (i.e. if we have only one
    # event above threshold, it might be an outlier in the noise, and if we recompute the threshold for this data i
    # it will find a threshold in the middle of background. this would estimate atom loading near 50%, even if it
    # was really 0%).
    over_thresh_fraction = sum([x > self.atom_counts_threshold for x in data])/len(data)

    if over_thresh_fraction > 0.3:
        # recompute the atoms_loaded using otsu thresholding. the optimizer will sometimes realize that it can
        # improve the cost by simply lowering the beam power until the threshold cuts the noise of an atom signal,
        # such that one atom trapped passes the threshold both ways several times before it leaves. This will be
        # incorrectly counted as many atoms, so to guard against this, we recompute the threshold from the data
        # using the otsu method.

        # logging.debug(f'initially computed {atoms_loaded} atoms loaded')
        threshold = threshold_otsu(data)
        q_last = (data[0] > threshold)
        for x in data[1:]:
            q = x > threshold
            if q != q_last and q_last:
                atoms_loaded += 1
            q_last = q
        atoms_loaded += q_last
        # logging.debug(f'recomputed computed {atoms_loaded} atoms loaded')

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

    shot1 = self.counts_list
    atoms_loaded = [x > self.single_atom_counts_threshold for x in shot1]
    n_atoms_loaded = sum(atoms_loaded)
    loading_fraction = n_atoms_loaded/len(shot1)
    return -100 * loading_fraction


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

    shot1 = self.counts_list
    atoms_loaded = [x > self.single_atom_counts_threshold for x in shot1]
    n_atoms_loaded = sum(atoms_loaded)
    loading_fraction = n_atoms_loaded/len(shot1)

    # If there are enough atoms loaded according to the threshold, recompute the loading rate with an Otsu threshold.
    # this will typically give a more accurate cut-off in case the histogram cleanness or cut-off changes with the
    # parameters we are varying. The reason we can not start with Otsu thresholding right away is that it would
    # still return a cut-off even if we load no atoms, and the cut-off would just bisect the background mode. Put
    # another way, it can not tell whether the data is bimodal or not.
    if loading_fraction > 0.3:  # apparent very low rate loading might just be wrongly classified background
        threshold = threshold_otsu(np.array(self.counts_list))
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
    shot1 = self.counts_list
    shot2 = self.counts2_list
    atoms_loaded = [x > self.single_atom_counts_threshold for x in shot1]
    n_atoms_loaded = sum(atoms_loaded)
    atoms_retained = [x > self.single_atom_counts2_threshold and y for x, y in zip(shot2, atoms_loaded)]
    retention_fraction = 0 if not n_atoms_loaded > 0 else sum(atoms_retained) / n_atoms_loaded
    loading_fraction = n_atoms_loaded/len(shot1)

    if loading_fraction > 0.3:  # apparent very low rate loading might just be wrongly classified background
        threshold = threshold_otsu(np.array(self.counts_list))
        atoms_loaded = [x > threshold for x in shot1]
        n_atoms_loaded = sum(atoms_loaded)
        loading_fraction = n_atoms_loaded / len(shot1)
        retention_fraction = 0 if not n_atoms_loaded > 0 else sum(atoms_retained) / n_atoms_loaded
        cost *= threshold/500

    # 0.6 is probably the best loading rate we can hope for.
    cost *= -100 * retention_fraction * loading_fraction / 0.6
    return cost


def atom_retention_cost(self) -> TFloat:
    """
    the cost function for optimizing the fraction of atoms retained in a two-shot experiment.

    :param self: experiment instance
    :return: -100*retention_fraction, the negated percentage of atoms detected in the readout
    """

    shot1 = self.counts_list
    shot2 = self.counts2_list
    atoms_loaded = [x > self.single_atom_counts_threshold for x in shot1]
    n_atoms_loaded = sum(atoms_loaded)
    atoms_retained = [x > self.single_atom_counts2_threshold and y for x, y in zip(shot2, atoms_loaded)]
    retention_fraction = 0 if not n_atoms_loaded > 0 else sum(atoms_retained) / n_atoms_loaded
    loading_fraction = n_atoms_loaded/len(shot1)

    # recompute the retention and loading with an Otsu threshold.
    # this will typically give a more accurate cut-off in case the histogram cleanness or cut-off changes with the
    # parameters we are varying. The reason we can not start with Otsu thresholding right away is that it would
    # still return a cut-off even if we load no atoms, and the cut-off would just bisect the background mode. Put
    # another way, it can not tell whether the data is bimodal or not.
    if loading_fraction > 0.3:  # apparent very low rate loading might just be wrongly classified background
        threshold = threshold_otsu(np.array(self.counts_list))
        atoms_loaded = [x > threshold for x in shot1]
        n_atoms_loaded = sum(atoms_loaded)
        atoms_retained = [x > self.single_atom_counts2_threshold and y for x, y in zip(shot2, atoms_loaded)]
        retention_fraction = 0 if not n_atoms_loaded > 0 else sum(atoms_retained) / n_atoms_loaded

    return -100 * retention_fraction


def atom_blowaway_cost(self) -> TFloat:
    """
    the cost function for minimizing the fraction of atoms retained in a two-shot experiment.

    :param self: experiment instance
    :return: 100*(retention_fraction-1)
    """

    shot1 = self.counts_list
    shot2 = self.counts2_list
    atoms_loaded = [x > self.single_atom_counts_threshold for x in shot1]
    n_atoms_loaded = sum(atoms_loaded)
    atoms_retained = [x > self.single_atom_counts2_threshold and y for x, y in zip(shot2, atoms_loaded)]
    retention_fraction = 0 if not n_atoms_loaded > 0 else sum(atoms_retained) / n_atoms_loaded

    return 100*(retention_fraction - 1)
