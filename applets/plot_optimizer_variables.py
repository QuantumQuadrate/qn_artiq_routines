"""
Plot the values for the variables being optimized

applet command:
python "C:\...\qn_artiq_routines\applets\plot_optimizer_variables.py"
optimizer_vars_dataset optimizer_var1 --var2 --var3 --var4
"""

import numpy as np
import PyQt5  # make sure pyqtgraph imports Qt5
from PyQt5.QtCore import QTimer
import pyqtgraph
import matplotlib.colors as mcolors
import numpy as np

from artiq.applets.simple import TitleApplet

MAX_VARIABLE_NUMBER = 10  # if you really want to optimize more values than this, go for it


def generate_colorblind_friendly_colors(n):
    # Generate equally spaced hues in HSL color space
    hues = np.linspace(0, 360, n + 1)[:-1]

    # Convert HSL to RGB
    rgb_colors = 255*np.array([mcolors.hsv_to_rgb((hue / 360, 0.7, 0.8)) for hue in hues])

    # return rgb_colors
    return [(r, g, b) for r, g, b in rgb_colors]

class XYPlot(pyqtgraph.PlotWidget):
    def __init__(self, args):
        pyqtgraph.PlotWidget.__init__(self)
        self.args = args

        self.symbols = ['o']*MAX_VARIABLE_NUMBER
        self.colors = generate_colorblind_friendly_colors(MAX_VARIABLE_NUMBER)

    def data_changed(self, data, mods, title):
        try:
            var_names = data[self.args.var_names][1]
            n_variables = min(len(var_names), MAX_VARIABLE_NUMBER)

            print(var_names)

            # should be a numpy array where each row is a different data channel
            optimizer_var_data = []
            for i in range(n_variables):
                optimizer_var_data.append(data[getattr(self.args, f"var{i}")][1][1:])

        except KeyError:
            return

        all_updated = True
        # check that all datasets are the same length
        for i in range(n_variables-1):
            if len(optimizer_var_data[i+1]) != len(optimizer_var_data[i]):
                all_updated = False
                break

        # only updating when all are updated makes managing x pts easier
        if all_updated:

            x = np.arange(len(optimizer_var_data[0]))

            self.clear()
            for i in range(n_variables):
                self.plot(x, optimizer_var_data[i],
                          pen=self.colors[i],
                          symbol=self.symbols[i],
                          symbolBrush=self.colors[i],
                          symbolPen='w',
                          name=var_names[i])
            self.setTitle(title)
            self.addLegend()


def main():
    applet = TitleApplet(XYPlot)
    applet.add_dataset("var_names", "a list of names of the variables being optimized")
    applet.add_dataset("var0", "optimizer variable 1")
    for i in range(1, MAX_VARIABLE_NUMBER):
        applet.add_dataset(f"var{i}", f"optimizer variable {i}", required=False)
    applet.add_dataset("var_bounds", "list of tuples specifying var bounds. "
                                     "if given, the data will be normalized to the bounds", required=False)

    applet.run()


if __name__ == "__main__":
    main()
