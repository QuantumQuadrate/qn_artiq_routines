"""
arguments are counts, counts2, and a defined threshold to discriminate between atom
and background.

applet command:
python "C:\..\qn_artiq_routines\applets\plot_retention_and_loading.py"
SPCM0_RO1 SPCM0_RO2 n_measurements iteration single_atom_threshold t_SPCM_first_shot

--passes_per_point <dataset> (optional): the shots each scan point holds, one
entry per point after a seed. Two_nodes_atom_atom_parity_experiment needs it:
its n_measurements counts heralds, but the shots are written once per pass.
Without it, or while it holds only its seed, the shots are grouped by
n_measurements.
"""

#!/usr/bin/env python3

import numpy as np
import PyQt5  # make sure pyqtgraph imports Qt5
from PyQt5.QtCore import QTimer
import pyqtgraph
from skimage.filters import threshold_otsu
from artiq.applets.simple import TitleApplet


def point_slices(n_shots, measurements, passes_per_point=None):
    """(start, stop) of each complete scan point's shots.

    passes_per_point: the shots of each point (seed already dropped), or None
    / empty to use measurements per point.
    """
    if passes_per_point is not None and len(passes_per_point):
        stops = np.cumsum(np.asarray(passes_per_point, dtype=int))
        stops = stops[stops <= n_shots]
        starts = np.concatenate([[0], stops[:-1]]).astype(int)
        return list(zip(starts, stops))
    return [(i*measurements, (i+1)*measurements) for i in range(n_shots//measurements)]


class XYPlot(pyqtgraph.PlotWidget):
    def __init__(self, args):
        pyqtgraph.PlotWidget.__init__(self)
        self.args = args
        self.timer = QTimer()
        self.timer.setSingleShot(True)
        # self.timer.timeout.connect(self.length_warning)
        self.mismatch = {'X values': False,
                         'Error bars': False,
                         'Fit values': False}

    def data_changed(self, data, mods, title):

        try: # not all of these are persistent
            counts_shot1 = np.array(data[self.args.counts_shot1][1][1:])
            counts_shot2 = np.array(data[self.args.counts_shot2][1][1:])
            measurements = data[self.args.measurements][1]
            threshold_cts_per_s = data[self.args.threshold_cts_per_s][1]
            t_exposure = data[self.args.t_exposure][1]
            cutoff = int(t_exposure*threshold_cts_per_s)

            passes_per_point = None
            if getattr(self.args, "passes_per_point", None):
                try:
                    passes_per_point = np.array(data[self.args.passes_per_point][1][1:], dtype=int)
                except (KeyError, TypeError, ValueError):
                    passes_per_point = None
            slices = point_slices(len(counts_shot1), measurements, passes_per_point)
            iteration = len(slices)
            if iteration > 0:
                if len(counts_shot1) == len(counts_shot2):

                    retention_array = np.zeros(iteration)
                    loading_rate_array = np.zeros(iteration)
                    n_atoms_loaded_array = np.zeros(iteration)

                    x = np.arange(iteration)

                    try:
                        nsteps = len(data.get(self.args.scan_sequence1, (False, None))[1])
                        scan_sequence1 = data[self.args.scan_sequence1][1]
                        if nsteps > 1 or scan_sequence1 != [0.0]:
                            x = np.array(scan_sequence1[:iteration])
                    except: # len will fail if sequence is None
                        pass

                    # if iteration > 0:
                    for i, (start, stop) in enumerate(slices):
                        shot1 = counts_shot1[start:stop]
                        shot2 = counts_shot2[start:stop]
                        atoms_loaded = [x > cutoff for x in shot1]
                        n_atoms_loaded = sum(atoms_loaded)
                        loading_fraction = n_atoms_loaded / (stop - start)

                        # could be helpful if the readout background is drifting, but needs to be tested for reliability
                        # if loading_fraction > 0.3:  # apparent very low rate loading might be wrongly classified background
                        #     cutoff = threshold_otsu(shot1)
                        #     atoms_loaded = [x > cutoff for x in shot1]
                        #     n_atoms_loaded = sum(atoms_loaded)
                        #     loading_fraction = n_atoms_loaded / measurements

                        n_atoms_loaded_array[i] = n_atoms_loaded
                        atoms_retained = [x > cutoff and y for x, y in zip(shot2, atoms_loaded)]
                        loading_fraction = n_atoms_loaded / (stop - start)
                        loading_rate_array[i] = loading_fraction
                        retention_fraction = 0 if not n_atoms_loaded > 0 else sum(atoms_retained) / n_atoms_loaded
                        retention_array[i] = retention_fraction

                    error_standard = np.array([1/np.sqrt(n) if n > 0 else 0 for n in n_atoms_loaded_array])
                    error_binomial = np.array([np.sqrt(n * p * (1 - p)) / n if n > 0 else 0 for n, p in zip(n_atoms_loaded_array, retention_array)])

                    error = error_standard

                    self.clear()
                    if len(x) == len(retention_array) and len(x) == len(loading_rate_array):
                        self.plot(x, retention_array,
                                  pen=None,
                                  symbol='o',
                                  symbolBrush=(255, 0, 0),
                                  symbolPen='w',
                                  name='retention')
                        self.plot(x, loading_rate_array,
                                  pen=None,
                                  symbol='o',
                                  symbolBrush=(0, 100, 100),
                                  symbolPen='w',
                                  name='loading')

                        self.setYRange(-0.0, 1.0, padding=0)

                        title = str(data[self.args.scan_vars][1])
                        self.setTitle(title)

                        if error is not None:
                            # See https://github.com/pyqtgraph/pyqtgraph/issues/211
                            if hasattr(error, "__len__") and not isinstance(error, np.ndarray):
                                error = np.array(error)
                        errbars = pyqtgraph.ErrorBarItem(
                            x=x, y=retention_array, height=2*error, pen=(255, 0, 0)) # error should be +/- the std, hence 2*
                        self.addItem(errbars)
                        self.addLegend()
            else:
                self.clear()
        except:
            self.clear()

def main():
    applet = TitleApplet(XYPlot)
    applet.add_dataset("counts_shot1", "first readout photocounts")
    applet.add_dataset("counts_shot2", "second readout photocounts")
    applet.add_dataset("measurements", "how many measurements per data point")
    applet.add_dataset("threshold_cts_per_s", "the atom signal counts/s threshold")
    applet.add_dataset("t_exposure", "the atom readout exposure time")
    applet.add_dataset("scan_vars", "the names of the scan variable(s)", required=False)
    applet.add_dataset("scan_sequence1", "the scan steps", required=False)
    applet.add_dataset("passes_per_point", "shots per scan point, if not n_measurements", required=False)

    applet.run()

if __name__ == "__main__":
    main()