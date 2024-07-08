"""
Connect to a K10CR1 device, query it with an RPC call from the kernel, and check the latency
"""

from time import sleep, time
from pylablib.devices import Thorlabs  # for Kinesis instrument control
import nidaqmx as daq # todo: remove
from artiq.experiment import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelmax
import logging
import os, sys

class KinesisMotorWrapper:

    def __init__(self, *args, **kwargs):
        self.motor = Thorlabs.KinesisMotor(*args, **kwargs)

    @rpc
    def get_position(self) -> TFloat:
        sleep(1)
        position = self.motor.get_position()
        return position

    @rpc
    def get_position_query_time(self) -> TFloat:
        _ = self.motor.get_position()
        return time()

@rpc
def get_host_time() -> TFloat:
    return time()


class RotatorCommTest(EnvExperiment):

    def build(self):
        self.setattr_argument('K10CR1_SN', StringValue('55000759'))
        self.setattr_device('core')

    def prepare(self):

        sn = int(self.K10CR1_SN)

        try:
            # self.rotor = Thorlabs.KinesisMotor(conn=sn, scale='K10CR1')
            self.rotor = KinesisMotorWrapper(conn=sn, scale='K10CR1')
            print(f"opened K10CR1 {sn}")
        except Exception as e:
            logging.error(f"Failed to connect to K10CR1 {sn}! "
                          f"\n Check that you correctly typed the SN and that it is connected.")
            print(sn, e)
            raise

        self.rotor_degrees = 0.0


    @kernel
    def run(self):

        self.core.reset()

        start = now_mu()
        at_mu(start)

        # the clock seems to idle when this is called, so that the timeline does not progress
        self.rotor_degrees = self.rotor.get_position()
        end = now_mu()
        print("delta mu time (s):", (end-start)*1e-9)

        # how much time does an RPC take?
        t_avg = 0.0
        for i in range(1000):
            t1 = get_host_time()
            t2 = get_host_time()
            t_avg += t2 - t1
        t_avg /= 1000
        print("avg RPC time (s):", t_avg)

        # how much time does querying the Thorlabs motor take?
        t_avg = 0.0
        for i in range(1000):
            t1 = self.rotor.get_position_query_time()
            t2 = self.rotor.get_position_query_time()
            t_avg += t2 - t1
        t_avg /= 1000
        print("avg RPC K10CR1 query time (s):", t_avg)




