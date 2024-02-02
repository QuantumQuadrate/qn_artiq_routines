
from time import sleep, time
from pylablib.devices import Thorlabs # for Kinesis instrument control
import nidaqmx as daq
import nidaqmx.constants as daq_constants
from nidaqmx.errors import DaqError, DaqWarning
from nidaqmx.error_codes import DAQmxErrors, DAQmxWarnings
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from K10CR1.k10cr1 import K10CR1
# a dictionary specifying channels we want to optimize

class RotatorFeedbackChannel():

    def __init__(self,ch_name, dds_channel=1, rotator_sn=['55105674', '55000741'], dry_run=True, max_runs=10, leave_laser_on=False):
        self.dry_run = dry_run
        if self.dry_run:
            self.measure = self.measure_dryrun  # use this for testing since you don't have access to artiq hardware yet
        else:
            self.measure = self._measure

        self.detector_volts = 0.0  # what we want to maximize

        self.ser0 = rotator_sn[0]
        self.ser1 = rotator_sn[1]

        self.scl = True  # whether to use physical units - this apparently has no effect
        self.stage = [Thorlabs.KinesisMotor(conn = self.ser0, scale='K10CR1'),
                      Thorlabs.KinesisMotor(conn = self.ser1, scale='K10CR1')]



        if self.stage[0].is_homed() is False:
            self.stage[0].home()
        if self.stage[1].is_homed() is False:
            self.stage[1].home()

    def print_pos(self, rotor_num = 0):
        print("position = ", self._pos(rotor_num))

    def print_abs_pos(self, rotor_num = 0):
        print("position = ", self._abs_pos(rotor_num))

    def _pos(self, rotor_num = 0):
        return self.stage[rotor_num].get_position(scale=self.scl)
    def _abs_pos(self, rotor_num = 0):
        return self._pos(rotor_num) % 360

    #@kernel
    def _measure(self):
        return 0

    # measure and update self.detector_volts

    def measure_dryrun(self, Ip = 1, phase_psi = 45, phase_chi = 45 ):
        intensity = Ip * np.sin((np.pi * np.float64(self._pos(0) - phase_psi)) / (180))*\
                    np.cos((np.pi * np.float64(self._pos(1) - phase_chi)) / (180))
        return (intensity + 0.05*(np.random.rand()-0.5))
    # measure and update self.detector_volts

    def move_by(self, degrees, rotor_num = -1, velocity = None, r1 = None):
        r = r1
        if rotor_num == -1 and r is None:
            r = self.stage[0]
        else:
            r = self.stage[rotor_num]
        if velocity is None:
            r.move_by(degrees)
            r.wait_move()
        else:
            print(r.get_velocity_parameters())
            r.setup_velocity(max_velocity=velocity, scale=self.scl)
            print(r.get_velocity_parameters())
            r.move_by(degrees)
            r.wait_move()
        return 0
    def move_to(self, degrees, rotor_num = 0, velocity = None, r1 = None):
        r = r1
        if rotor_num == -1 and r is None:
            r = self.stage[0]
        else:
            r = self.stage[rotor_num]

        r.move_to(degrees)
        r.wait_move()
        return 1
    def jog(self, direction = "+", rotor_num = 0, velocity = 10):
        return 0
    def optimize(self, initial_velocity = 10, initial_acceleration = 10, starting_point = 0):
        """
        Set up initial velocity and acceleration parameters to be at their maximums
        in deg/s and deg/s^2
        """
        self.stage[0].setup_velocity(max_velocity=initial_velocity, acceleration=initial_acceleration,
                                  scale = self.scl)
        self.stage[1].setup_velocity(max_velocity=initial_velocity, acceleration=initial_acceleration,
                                     scale=self.scl)
        r0 = self.stage[0]
        r1 = self.stage[1]
        self.print_pos(0)
        self.print_pos(1)

        r0.move_to(0)
        r1.move_to(0)
        for i in range(0,2):
            print(f"Moving to {starting_point} degrees, and then waiting")
            self.move_to(starting_point, rotor_num = i)


        r0.wait_for_stop()
        r1.move_to(0)
        self.print_pos(0)
        self.print_pos(1)
        r0.move_to(0)
        r1.move_to(0)

        self.optimise_for_rotor(0, pos = 137, restrict_to_near=True, range = 10)
        self.optimise_for_rotor(1, pos = 42)
        self.close()
        # self.dds_channel.sw.on() # switch on the laser. don't worry about this
        """
        # the main loop
        for i in range(self.max_runs):
            self.measure()
            # based on self.detector_volts,
            self.move_rotator(1)
            changing = True
            tolerance_met = False
            result = False
            if tolerance_met or result is not changing:
                break"""

        # self.dds_channel.sw.off() # switch on the laser. don't worry about this
        return 0



    """
        # this would happen in BaseEperiment:
        class RotatorExperiment(EnvExperiment):
     
        def run(self):   
            for key, val in self.rotator_lib: #wrong
                # set a RotatorFeedbackChannel as attribute of our current experiment
                setattr(self.experiment, key + "rotator_ch", RotatorFeedbackChannel(**val, max_runs=10))
        def __init__(self):
            self.rotator_lib = {}
    """
    def close(self):
        self.stage[0].close()
        self.stage[1].close()
        print("Closed successfully")

    def optimise_for_rotor(self, rotor_num = 0, init_pos = 0, restrict_to_near = True, range = 45, favor_start = 1):
        stage = self.stage[rotor_num]
        stage.stop()
        stage.move_to(position=init_pos)
        stage.wait_for_stop()

        starting_time = time()

        stage.setup_jog(max_velocity=10, acceleration=10)
        stage.jog(direction="+")

        i = 0

        degree = stage.get_position()
        print(f"Starting position is {self._abs_pos(rotor_num=rotor_num)} degrees")
        # User can pass a float variable through favor_start so it takes a larger percentage change in intensity
        # to change it from it's original position

        max_degree = degree
        i_max = self.measure() * favor_start
        restriction_factor = 2 * favor_start



        if restrict_to_near:
            min_possible = init_pos - range
            max_possible = init_pos + range
        else:
            min_possible = 361
            max_possible = -361

        reversals = 0
        moving_forward = True

        time_wait = 0.1
        speed = 10

        while stage.is_moving():

            degree = stage.get_position()
            intensity = self.measure()

            if intensity >= i_max:
                i_max = intensity
                max_degree = degree
                was_above = True
                reversals = 0

            elif intensity <= i_max / restriction_factor or (degree >= max_possible or degree <= min_possible):
                restriction_factor /= 1.3

                if restriction_factor <= 1:
                    restriction_factor = 1.001

                stage.stop()

                if moving_forward is True:
                    stage.setup_jog(max_velocity=speed/(reversals+2), acceleration=10)
                    stage.jog(direction="-")
                    moving_forward = False
                    # max_possible = degree
                    sleep(2)
                    """while d > degree - (45.0 / 2.0 ** reversals):
                        d = stage.get_position()
                        sleep(0.1)"""
                else:
                    stage.setup_jog(max_velocity=speed / (reversals + 2), acceleration=10)
                    stage.jog(direction="+")
                    moving_forward = True
                    # min_possible = degree
                    sleep(2)
                    """while d < (degree + 45.0 / 2.0 ** reversals):
                        d = stage.get_position()
                        sleep(0.1)"""

                reversals += 1
                time_wait = (0.2 / (reversals))

                if reversals == 2:
                    stage.move_to(position=max_degree)
                    stage.wait_for_stop()

            sleep(time_wait)
            i += 1
            if i == 1000:
                stage.stop()
        # If the user passed a factor to favor the starting position it is factored out
        if favor_start != 1:
            i_max/=favor_start

        print(f"Rotor number{rotor_num}\n--------------------------")
        print(f"Max degree:{max_degree}\nIntensity at degree:{i_max}")
        print(f"time to optimize rotor number {rotor_num}: ({time()-starting_time} seconds)")
        stage.wait_for_stop()

def test():
    rotator_feedback_dict = {
        'ch_1':
            {
                "ch_name": "dipole_trap1",
                "dds_channel": "1",  # this is basically the laser we want to turn on. don't worry about this
                "rotator_sn": 0o012345,  # the rotator id or serial number
                "sampler_ch": "1"  # the channel which reads in the voltage from a detector. don't worry about this
            },
        'ch_2':
            {
                "ch_name": "dipole_trap2",
                "dds_channel": "2",
                "rotator_sn": 0o023456,  # the rotator id or serial number
                "sampler_ch": "2"  # the channel which reads in the voltage from a detector. don't worry about this
            }
    }
    def intensity_sim(degree, phase, A = 1):
        intensity = A*np.sin((np.pi*np.float64(degree-phase))/(180))
        return intensity

    """"""
    stage = Thorlabs.KinesisMotor(conn="55105674", scale='K10CR1')
    #stage.move_by(100)
    #stage.wait_move()
    stage.stop()
    stage.wait_for_stop()
    sleep(1)
    phase = 45

    stage.move_to(position=132)
    stage.wait_for_stop()
    stage.jog(direction = "+")

    i = 0

    degree = stage.get_position()
    print(degree)
    i_max = intensity_sim(degree=degree, phase=phase)
    max_degree = -361
    min_possible = -361
    max_possible = 180
    reversals = 0
    moving_forward = True
    factor = 1.5
    time_wait = 0.1
    while stage.is_moving():
        degree = stage.get_position()
        intensity = intensity_sim(degree=degree, phase=phase)

        if intensity >= i_max:
            i_max = intensity
            max_degree = degree
            was_above = True

        elif intensity <= i_max/factor or (degree >= max_possible or degree <= min_possible):
            factor/=1.1
            if factor <= 1:
                factor = 1.001

            stage.stop()

            if moving_forward is True:
                stage.jog(direction="-")
                moving_forward = False
                #max_possible = degree
                sleep(2)
                """while d > degree - (45.0 / 2.0 ** reversals):
                    d = stage.get_position()
                    sleep(0.1)"""
            else:
                stage.jog(direction="+")
                moving_forward = True
                #min_possible = degree
                sleep(2)
                """while d < (degree + 45.0 / 2.0 ** reversals):
                    d = stage.get_position()
                    sleep(0.1)"""

            reversals += 1
            time_wait = (0.2/(reversals))
            if reversals == 10:

                stage.move_to(position=max_degree)
                stage.wait_for_stop()

        sleep(time_wait)
        i+=1
        if i == 10000:
            stage.stop()

    print(f"Max degree:{max_degree}\nIntensity at degree:{i_max}")
    stage.close()

      # don't do anything until the motion is complete
    #rotor1.print_pos()
    #rotor1.test()

def main():


    devices = Thorlabs.list_kinesis_devices()  # can use this to print out the detected devices
    print(devices)
    rotor1 = RotatorFeedbackChannel(ch_name="dipole_trap1", rotator_sn=["55105674", "55000741"])
    """print(rotor1.stage[0].get_all_channels())
    rotor1.print_pos()
    ref_degs = 0  # the angle of the polarizer about which to move +/-
    rotor1.print_pos()
    print(rotor1.move_by(10,))
    sleep(10)
    rotor1.move_by(-10)"""
    rotor1.stage[0].setup_velocity(max_velocity=10)
    rotor1.optimize()
    rotor1.close()

if __name__ == '__main__':
    main()






