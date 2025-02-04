"""
For plotting up to 10 numeric sequences, for example, the DDS RF power history

Example usage:
python "C:\..\qn_artiq_routines\applets\plot_xy_multichannel.py"
p_AOM_A1_history MOT_beam_monitor_points --y2 p_AOM_A2_history --y3 p_AOM_A3_history
--y4 p_AOM_A4_history --y5 p_AOM_A5_history --y6 p_AOM_A6_history --y7 p_FORT_loading_history --labels feedbackchannels
"""

import numpy as np
import PyQt5  # make sure pyqtgraph imports Qt5
from PyQt5.QtCore import QTimer
import pyqtgraph
from matplotlib import pyplot as plt

from artiq.applets.simple import TitleApplet


class XYPlot(pyqtgraph.PlotWidget):
    def __init__(self, args):
        pyqtgraph.PlotWidget.__init__(self)
        self.args = args
        self.timer = QTimer()
        self.timer.setSingleShot(True)
        self.timer.timeout.connect(self.length_warning)
        self.mismatch = {'X values': False,
                         'Error bars': False,
                         'Fit values': False}

    def average_over_measurement(self, measurements, iteration, history):
        mean_by_iteration = [np.mean(history[j * measurements:(j + 1) * measurements]) for j in range(iteration)]
        return mean_by_iteration


    def data_changed(self, data, mods, title):
        try:
            measurements = data[self.args.measurements][1]

            y1 = np.array(data[self.args.y1][1][1:])
            iteration = len(y1) // measurements

            labels = data.get(self.args.labels, (False, None))[1]

            colors = plt.rcParams["axes.prop_cycle"]()

            if labels is None:
                labels = list(range(10))

            self.clear()
            for i in range(10):
                try:
                    y_data = data[getattr(self.args, f"y{i+1}")][1]
                    mean_by_iteration = self.average_over_measurement(measurements, iteration, y_data)
                    self.plot(range(iteration), mean_by_iteration, pen=next(colors)["color"], name=labels[i])
                except KeyError as e:
                    break

            self.addLegend()

        except KeyError:
            return

    def length_warning(self):
        self.clear()
        text = "⚠️ dataset lengths mismatch:\n"
        errors = ', '.join([k for k, v in self.mismatch.items() if v])
        text = ' '.join([errors, "should have the same length as Y values"])
        self.addItem(pyqtgraph.TextItem(text))


def main():
    applet = TitleApplet(XYPlot)
    applet.add_dataset("y1", "y1 values")
    applet.add_dataset("measurements", "how many measurements per data point")
    applet.add_dataset("y2", "y2 values", required=False)
    applet.add_dataset("y3", "y3 values", required=False)
    applet.add_dataset("y4", "y4 values", required=False)
    applet.add_dataset("y5", "y5 values", required=False)
    applet.add_dataset("y6", "y6 values", required=False)
    applet.add_dataset("y7", "y7 values", required=False)
    applet.add_dataset("y8", "y8 values", required=False)
    applet.add_dataset("y9", "y9 values", required=False)
    applet.add_dataset("y10", "y10 values", required=False)
    applet.add_dataset("x", "x values", required=False)
    applet.add_dataset("labels", "channel names", required=False)
    # applet.add_dataset("error", "Error bars for each X value", required=False)
    # applet.add_dataset("fit", "Fit values for each X value", required=False)
    applet.run()


if __name__ == "__main__":
    main()
