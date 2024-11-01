from time import time, sleep
from pylablib.devices import Thorlabs  # for Kinesis instrument control
from artiq.experiment import *

# driver class. We'll expose the functions of this class
# after it's been instantiated



class K10CR1_NDSP_Driver():

    def __init__(self, **kwargs):
        self.motor = Thorlabs.KinesisMotor(**kwargs)

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
    
class K10CR1_Multi_NDSP_Driver():
    """
    Wrapper class to allow multiple K10CR1 units to be controlled
    with one NDSP server.

    I want to be able to do something like:
    self.k10cr1_ndsp.move_by(["780_QWP","780_HWP"],
                             [10, 20])

    which implies the move_by method operates like this:
    def move_by(self, names, degrees):
        for name,deg in zip(names, degrees):
            self.rotors[name].move_by(deg)

    can I generalize this? don't try to be too clever
    """
    def __init__(self, devices_kwargs):
        """
        devices_kwargs is a list of tuples of strings and dictionaries:
        [(name1, kwargs1), (name1, kwargs2), ...]
        where kwargs1, etc will be used to set up the device instance
        """

        self.motor_list = {}

        for name, kwargs in devices_kwargs:
            self.motor_list[name] = Thorlabs.KinesisMotor(**kwargs)

    @rpc
    def get_position(self, names) -> TArray(TFloat):
        positions = []
        for name in names:
            assert name in self.motor_list, f"{name} does not exist! did you mistype it?"
            positions.append(self.motor_list[name].get_position())
        return positions

    @rpc
    def is_moving(self, names) -> TBool:
        """
        returns true if any motors are moving
        """
        moving = 0
        for name in names:
            assert name in self.motor_list, f"{name} does not exist! did you mistype it?"
            if self.motor_list[name].is_moving():
                moving = 1
                print(f"K10CR1 {name} moving!")
        return moving

    @rpc(flags={'async'})
    def move_by(self, names, degrees):
        for name,deg in zip(names, degrees):
            assert name in self.motor_list, f"{name} does not exist! did you mistype it?"
            self.motor_list[name].move_by(deg)
            print(f"K10CR1 {name} moving!")

    @rpc(flags={'async'})
    def move_to(self, names, degrees):
        for name,deg in zip(names, degrees):
            assert name in self.motor_list, f"{name} does not exist! did you mistype it?"
            self.motor_list[name].move_by(deg)
            print(f"K10CR1 {name} moving!")

    @rpc(flags={'async'})
    def wait_move(self, names):
        for name in names:
            assert name in self.motor_list, f"{name} does not exist! did you mistype it?"
            self.motor_list[name].wait_move()
            print(f"K10CR1 {name} moving!")

