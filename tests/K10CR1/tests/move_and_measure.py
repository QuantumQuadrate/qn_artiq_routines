from artiq.experiment import *
import numpy as np
from pylablib.devices import Thorlabs  # for Kinesis instrument control
import matplotlib.pyplot as plt

import sys, os
cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")

from K10CR1.KinesisMotorWrapper import KinesisMotorWrapper


class MoveK10CR1andMeasure(EnvExperiment):

    def build(self):
        self.setattr_device("core")
        self.setattr_device("sampler1")

    def prepare(self):

        self.sampler_buffer = np.zeros(8)
        self.sampler_ch = 6  # channel to measure

        self.rotor_SN = self.get_dataset("K10CR1_852_HWP_SN")

    @rpc
    def move_and_wait(self):
        # the rotor instance is destroyed at the end of each rpc call, so make it fresh each time
        # rotor = Thorlabs.KinesisMotor(conn=, **kwargs)
        pass