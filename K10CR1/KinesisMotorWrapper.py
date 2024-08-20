from time import time, sleep
from pylablib.devices import Thorlabs  # for Kinesis instrument control
from artiq.experiment import *
import numpy as np


class KinesisMotorWrapper:
    """
    wrapper class to explicitly declare return types and add some functions for diagnostics

    Where applicable, the methods here shadow the names of the methods belonging to Thorlabs.KinesisMotor

    For now, the hardware settings can be configured in the APT User program with 'persist to hardware' checked.
    See ./setup
    """

    def __init__(self, *args, **kwargs):
        self.motor = Thorlabs.KinesisMotor(*args, **kwargs)

    @rpc
    def get_position(self) -> TFloat:
        position = self.motor.get_position()
        return position

    @rpc
    def is_moving(self) -> TBool:
        moving = self.motor.is_moving()
        return moving

    @rpc(flags={'async'})
    def move_by(self, degrees):
        self.motor.move_by(degrees)
        print("K10CR1 moving!")

    @rpc(flags={'async'})
    def move_to(self, degrees):
        self.motor.move_to(degrees)

    @rpc(flags={'async'})
    def wait_move(self):
        self.motor.wait_move()

    # additional functions for diagnostics
    @rpc
    def get_position_query_time(self) -> TFloat:
        _ = self.motor.get_position()
        return time()


class KinesisMotorSimulator:
    """
    A simulator of the KinesisMotorWrapper class for hardware-free tests

    Note that this is a rather simple class and does not attempt to mimic the full the thing
    except in very basic ways. Methods such as is_moving, e.g., will always return False as
    I haven't taken the time to implement threading so we can use non-blocking function calls.

    Where applicable, the methods here shadow the names of the methods belonging to Thorlabs.KinesisMotor
    """

    def __init__(self, *args, **kwargs):
        self.position = np.pi

    @rpc
    def get_position(self) -> TFloat:
        sleep(0.01)
        return self.position

    @rpc
    def is_moving(self) -> TBool:
        sleep(0.01)
        return False

    @rpc(flags={'async'})
    def move_by(self, degrees):
        sleep(0.1)
        self.position += degrees

    @rpc(flags={'async'})
    def move_to(self, degrees):
        sleep(.1)
        self.position = degrees

    @rpc(flags={'async'})
    def wait_move(self):
        sleep(0.1)

    # additional functions for diagnostics
    @rpc
    def get_position_query_time(self) -> TFloat:
        sleep(0.01)
        return time()
