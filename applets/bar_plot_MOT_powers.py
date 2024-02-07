"""
Plot the signals representing the laser beam powers. units are determined by the datasets

applet command:
python "C:\...\qn_artiq_routines\applets\bar_plot_MOT_powers.py"
MOT1_monitor MOT2_monitor MOT3_monitor MOT4_monitor MOT5_monitor MOT6_monitor
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
        self.labels = [f'MOT{i + 1}' for i in range(6)]
        # red-green color-blind friendly RGB colors from ChatGPT
        self.colors = [(0, 92, 169),  # blue
                        (255, 138, 0),  # orange
                        (149, 0, 210),  # purple
                        (0, 170, 150),  # teal
                        (128, 128, 128),  # gray
                        (139, 69, 19),  # brown
                        (255, 105, 180)]  # pink


    def data_changed(self, data, mods, title):
        try:
            # the data display will be rolling, only showing display_pts at a time

            # should be a numpy array where each row is a different data channel
            MOT_data = []
            for i in range(6):
                MOT_data.append(data[getattr(self.args, self.labels[i])][1][-1])

            MOT_switchyard_input = data.get(self.args.MOT_switchyard_input, (False, None))[1]

            if MOT_switchyard_input is not None:
                MOT_data.append(data[self.args.MOT_switchyard_input][1][-1])
                self.labels += ['MOT_switchyard_input']

        except KeyError:
            return

        x = range(len(MOT_data))

        self.clear()
        w = 0.8
        xpts = range(-1,len(MOT_data)+1)
        bar_graph = pyqtgraph.BarGraphItem(x = x, height = MOT_data, width = w, brushes=self.colors)
        self.setYRange(0.0, 1.1, padding=0)
        self.setXRange(-w, len(MOT_data)-1+w, padding=0)
        self.addItem(bar_graph)
        self.plot(xpts, np.full(len(xpts), 1), pen='red', style=PyQt5.QtCore.Qt.DashLine, width=0.2)
        self.setTitle(title)

        for i,datum in enumerate(MOT_data):
            if datum > 0.5:
                self.text = pyqtgraph.TextItem(str(round(datum,3)),color=(0,0,0))
                self.addItem(self.text)
                self.text.setPos(i-w/3, round(datum,3)-0.05)
            else:
                self.text = pyqtgraph.TextItem(str(round(datum, 3)), color=(255, 255, 255))
                self.addItem(self.text)
                self.text.setPos(i - w / 3, round(datum, 3) + 0.1)

def main():
    applet = TitleApplet(XYPlot)
    applet.add_dataset("MOT1", "MOT1 PD voltage")
    applet.add_dataset("MOT2", "MOT2 PD voltage")
    applet.add_dataset("MOT3", "MOT3 PD voltage")
    applet.add_dataset("MOT4", "MOT4 PD voltage")
    applet.add_dataset("MOT5", "MOT5 PD voltage")
    applet.add_dataset("MOT6", "MOT6 PD voltage")
    applet.add_dataset("MOT_switchyard_input", "MOT PD0 voltage", required=False)
    applet.run()

if __name__ == "__main__":
    main()
