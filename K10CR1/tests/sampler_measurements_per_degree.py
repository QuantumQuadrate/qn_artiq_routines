from artiq.experiment import *
import logging
import numpy as np
from scipy.optimize import curve_fit, minimize
import matplotlib.pyplot as plt
from time import time

import sys, os
cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")

from subroutines.rotator_feedback import FORTPolarizationOptimizer


class SamplerMeasurementsPerDegree(EnvExperiment):

    def build(self):
        self.setattr_device("core")
        self.setattr_device("sampler1")
        # self.setattr_argument("is_this_working", BooleanValue(False))

    def prepare(self):

        self.sampler_buffer = np.zeros(8)
        self.sampler_ch = 6  # channel to measure

        self.n = 10
        self.degrees_to_move = np.linspace(1, 90, self.n)
        self.measurements_list = [0] * self.n

        self.hwp_measurements = [0] * self.n
        self.qwp_measurements = [0] * self.n

        self.pol_optimizer = FORTPolarizationOptimizer(
            experiment=self, sampler=self.sampler1, sampler_ch=6, max_moves=10, HWP_SN='55000759', QWP_SN='55000740',
            tolerance=0.05, debugging=True, dry_run=False
        )

    @kernel
    def get_hwp_measurements(self) -> TList(TInt32):
        self.core.reset()

        for i in range(self.n):
            self.measurements_list[i] = self.pol_optimizer.sampler_pts_per_angle(self.degrees_to_move[i], 'hwp')

        return self.measurements_list

    @kernel
    def get_qwp_measurements(self) -> TList(TInt32):
        self.core.reset()

        for i in range(self.n):
            self.measurements_list[i] = self.pol_optimizer.sampler_pts_per_angle(self.degrees_to_move[i], 'qwp')

        return self.measurements_list

    def run(self):
        self.hwp_measurements = self.get_hwp_measurements()
        # self.qwp_measurements = self.get_qwp_measurements()
        print("done making measurements")

    def analyze(self):
        print("now we're plotting?")

        plt.plot(self.degrees_to_move, [h/deg for h, deg in zip(self.hwp_measurements, self.degrees_to_move)], label='HWP')
        # plt.plot(range(self.n), [q/deg for q, deg in zip(self.qwp_measurements, self.degrees_to_move)], label='HWP')
        plt.ylabel("measurements per\ndegrees moved")
        plt.show()