"""
Plot the signals representing the laser beam powers. units are determined by the datasets

applet command:
python "C:\...\qn_artiq_routines\applets\plot_MOT_powers.py"
MOT1_monitor MOT2_monitor MOT3_monitor MOT4_monitor MOT5_monitor MOT6_monitor
MOT_switchyard_monitor MOT_beam_monitor_points
"""

import numpy as np
import PyQt5  # make sure pyqtgraph imports Qt5
from PyQt5.QtCore import QTimer
import pyqtgraph

from artiq.applets.simple import TitleApplet


class XYPlot(pyqtgraph.PlotWidget):
    def __init__(self, args):
        pyqtgraph.PlotWidget.__init__(self)
        self.args = args
        self.labels = [f'MOT{i + 1}' for i in range(6)] + ['MOT_switchyard_input']
        # red-green color-blind friendly RGB colors from ChatGPT
        self.colors = [(0, 92, 169),  # blue
                        (255, 138, 0),  # orange
                        (149, 0, 210),  # purple
                        (0, 170, 150),  # teal
                        (128, 128, 128),  # gray
                        (139, 69, 19),  # brown
                        (255, 105, 180)]  # pink

        self.symbols = ['o', 't', 's', 't2', 'h', 't1', 'd']

    def data_changed(self, data, mods, title):
        try:
            # the data display will be rolling, only showing display_pts at a time
            pts = (data[self.args.pts][1])

            # should be a numpy array where each row is a different data channel
            MOT_data = []
            for i in range(6):
                MOT_data.append(data[getattr(self.args, self.labels[i])][1][-pts:])
            MOT_data.append(data[self.args.MOT_switchyard_input][1][-pts:])

        except KeyError:
            return

        all_updated = True
        # check that all datasets are the same length
        for i in range(6):
            if len(MOT_data[i+1]) != len(MOT_data[i]):
                all_updated = False
                break

        # only updating when all are updated makes managing x pts easier
        if all_updated:

            x = data.get(self.args.x, (False, None))[1]

            if x is None:
                x = np.arange(len(MOT_data[0]))
            else:
                x = x[-pts:]

            self.clear()
            for i in range(7):
                self.plot(x, MOT_data[i],
                          pen=self.colors[i],
                          symbol=self.symbols[i],
                          symbolBrush=self.colors[i],
                          symbolPen='w',
                          name=self.labels[i])
            self.setTitle(title)
            self.addLegend()
            # todo: use timestamps on the x axis?
            #  axis = DateAxisItem()
            #  plot.setAxisItems({'bottom':axis})

def main():
    applet = TitleApplet(XYPlot)
    applet.add_dataset("MOT1", "MOT1 fW PD voltage")
    applet.add_dataset("MOT2", "MOT2 fW PD voltage")
    applet.add_dataset("MOT3", "MOT3 fW PD voltage")
    applet.add_dataset("MOT4", "MOT4 fW PD voltage")
    applet.add_dataset("MOT5", "MOT5 fW PD voltage")
    applet.add_dataset("MOT6", "MOT6 fW PD voltage")
    applet.add_dataset("MOT_switchyard_input", "MOT PD0 voltage")
    applet.add_dataset("pts", "number of points to display")
    applet.add_dataset("x", "X values", required=False)
    applet.run()

if __name__ == "__main__":
    main()
