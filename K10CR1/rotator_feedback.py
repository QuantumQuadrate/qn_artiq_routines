
from time import sleep, time
from pylablib.devices import Thorlabs # for Kinesis instrument control
import nidaqmx as daq
from artiq.experiment import *
import nidaqmx.constants as daq_constants
from nidaqmx.errors import DaqError, DaqWarning
from nidaqmx.error_codes import DAQmxErrors, DAQmxWarnings
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelmax
import logging
# a dictionary specifying channels we want to optimize

def numerical_partial_derivative_x(function, x, y, epsilon=0.5):

    return (function(x + epsilon, y) - function(x - epsilon, y)) / (2 * epsilon)


def numerical_partial_derivative_y(function, x, y, epsilon=0.5):

    return (function(x, y + epsilon) - function(x, y - epsilon)) / (2 * epsilon)
class RotatorFeedbackChannel():

    def __init__(self,ch_name = "Dev1/ai0", dds_channel=0, rotator_sn=['55000741','55105674',], dry_run = True ,
                 max_runs=10, spc = None, rate = None,):

        self.spc = spc
        self.dry_run = dry_run
        self.daq_task = daq.Task()
        self.rate = rate

        if self.spc is None:
            self.spc = 1

        if self.rate is None or self.rate > 1e5:
            self.rate = 1e5


        self.daq_task.ai_channels.add_ai_voltage_chan(physical_channel=ch_name)

        if self.dry_run:
            self.measure = None
            print("No dry run implemented, edit code")
            # use this for testing since you don't have access to artiq hardware yet
        else:
            self.measure = self._measure

        self.ser0 = rotator_sn[0]
        self.ser1 = rotator_sn[1]

        self.scl = True  # whether to use physical units - this apparently has no effect
        self.stage = [Thorlabs.KinesisMotor(conn=self.ser0, scale='K10CR1', ),
                      Thorlabs.KinesisMotor(conn=self.ser1, scale='K10CR1', )]
        self.extra_rotors = []
        if (len(rotator_sn) > 2):
            for i in range(2, len(rotator_sn)):
                new_rotor = Thorlabs.KinesisMotor(conn=rotator_sn[i], scale = 'K10CR1')
                self.extra_rotors.append(new_rotor)
        self.r0 = self.stage[0]
        self.r1 = self.stage[1]
        self.r0.set_position_reference()
        self.r1.set_position_reference()

        #(f"{self._pos(0)} : {self._pos(1)}; intensity:{np.sum(self.daq_task.read(number_of_samples_per_channel=1))}")

    def print_connections(self, all = True, rotors = False, task = False):
        if rotors is True:
            for r in self.stage:
                print(r)
            for r in self.extra_rotors:
                print(r)

    def print_pos(self, rotor_num=0):
        print("position = ", self._pos(rotor_num))

    def wait_stop(self, rotor_num = 2):
        if rotor_num == 2:
            self.stage[0].wait_for_stop()
            self.stage[1].wait_for_stop()
        else:
            self.stage[rotor_num].wait_for_stop()

    def is_moving(self, rotor_num = None):
        if rotor_num is None:
            return self.r0.is_moving() or self.r1.is_moving()
        else:
            return self.stage[rotor_num].is_moving()

    def print_abs_pos(self, rotor_num=0):
        print("position = ", self._abs_pos(rotor_num))

    def _pos(self, rotor_num=0):
        return self.stage[rotor_num].get_position(scale=self.scl)

    def _abs_pos(self, rotor_num=0):
        return self._pos(rotor_num) % 180

    """Using the USB-NiDaq reads out voltage measured from a Si Amplified Detector"""
    def _measure(self, measurements=None):
        if measurements is None:
            measurements = self.spc

        data = self.daq_task.read(number_of_samples_per_channel=measurements)

        return np.mean(data)




    def q_h_gen(self, range, steps, center_x, center_y, random=False):
        """

        :param range: width and height of the range that the random/pre-generated points will be within
        :param steps: amount of points
        :param center_x: center of axis of the range for "rotor1"
        :param center_y: center of axis of the range for "rotor2"
        :param random: if False generates purely uniformly random points with specified range and center
        :return: x, y arrays of the random points where len(x) = len(y) = steps
        """
        if not random:
            range /= 2
            rads = np.arange(0, 2 * np.pi, (2 * np.pi) / steps)
            x = range * (np.cos(rads)) + center_x
            y = range * (np.sin(rads)) + center_y
            x1 = range / 2 * (np.cos(rads)) + center_x
            y1 = range / 2 * (np.sin(rads)) + center_y
            x2 = range / 4 * (np.cos(rads)) + center_x
            y2 = range / 4 * (np.sin(rads)) + center_y
            x0 = np.append(x, x1)
            x0 = np.append(x0, x2)
            y0 = np.append(y, y1)
            y0 = np.append(y0, y2)

        if random is True:
            x0 = range*(np.random.rand(steps) - 1/2) + center_x
            y0 = range*(np.random.rand(steps) - 1/2) + center_y

        return x0, y0


    def move_and_measure(self, range=None, steps=None, a=0, b=0,
                         dry_run=True, warnings = True):
        """
            Using pre selected point/randomly generated ones, the rotors are moved between points and then stopped.
            While the rotors are in the process of moving, we are keeping track of the angle it is at, and then the measurement.
            Due to time to measure, acceleration and speed of rotor, there is slight discrepancy between measured angle and
            actual angle for some given intensity measurement, but the non acceleration component of the discrepancy is solved
            for in the optimize portion of the code, and the acceleration error component seems to be negligible, when the goal
            is achieving maximum intensity.
            TODO sort points by distance before moving and measuring so that code will run quicker and more points can be included in path without
            changing time of computation
        """
        """
        :param: range: range of values in both directions that will be walked through
        :param: steps: the number of 2d points that will be moved between during the random walk
        :param: a = initial position of rotor1
        :param: b = initial position of rotor2
        :return: returns x, y, z of points where the first two are the postions and the third is the measured intensity  
        """
        ang1, ang2 = self.q_h_gen(steps=steps, range=range, center_x=a, center_y=b, random=True)
        m_ang1 = []
        m_ang2 = []
        measure = self.measure
        measurements = []
        if dry_run and warnings is True:
            ##TODO
            logging.warning("Dry run is unimplemented")
            return -1
        else:
            for a1, a2 in zip(ang1, ang2):
                self.r0.move_to(a1)
                self.r1.move_to(a2)
                i = 0
                sleep(0.3)
                while self.is_moving():
                    m_ang1 = np.append(m_ang1, self._pos(rotor_num=0))
                    m_ang2 = np.append(m_ang2, self._pos(rotor_num=1))
                    measurement = measure(measurements=1)
                    measurements = np.append(measurements, measurement)

                    i += 1
                    sleep(0.3)
                self.stage[1].wait_for_stop()
                self.stage[0].wait_for_stop()

        return m_ang1%180, m_ang2%180, measurements
    def move_by(self, degrees, rotor_num = -1, velocity = None, r1 = None):
        """
        Moves by
        :param degrees: distance to move can be +/- any float, not accurate below ~ 0.5 degs
        :param rotor_num: the number rotor to move
        :param velocity: speed to move the specified distance, max and default is 10 deg/s
        :param r1: any specified Thorlabs.KinesisMotor class variable
        :return:
        """
        r = r1
        if rotor_num == -1 and r is None:
            r = self.stage[0]
        elif r1 is None:
            r = self.stage[rotor_num]
        if velocity is None:
            r.move_by(degrees)
            r.wait_move()
        else:
            r.setup_velocity(max_velocity=velocity, scale=self.scl)
            r.move_by(degrees)
            r.wait_move()
        return 0

    def move_to(self, degrees, rotor_num, velocity = None, r1 = None):
        r = r1
        if rotor_num == -1 and r is None:
            logging.warning("Must provide rotor number or instantiated rotor class: no movement will be made")
            return -1
        else:
            r = self.stage[rotor_num]

        r.move_to(degrees)
        r.wait_move()
        return 1
    def function_to_maximize(self,x,y,):
        """

        :param x: the position that the first initialized rotor "self.stage[0]" will be moved to
        :param y: the position that the first initialized rotor "self.stage[0]" will be moved to
        :return: measurement value read for NiDaq
        """
        self.move_to(x, 0)
        self.move_to(y, 1)
        #ensures that the rotors are at the correct position before any measurement is made
        self.wait_stop(rotor_num=0)
        self.wait_stop(rotor_num=1)
        measurement = self.measure(measurements=1000)
        sleep(0.1) #sleeping to ensure that the measurement is completed before the rotors start moving again
        return measurement
    def gradient_ascent_2d(self, learning_rate, tolerance, min_iterations=10, max_iteration=100, init_x = 0, init_y = 0,
                           epsilon = 0.5):
        # Start with initial guesses for the maximum (x, y)
        x = init_x
        y = init_y
        z = self.function_to_maximize(x, y)
        xs = []
        ys = []
        zs = []
        new_gradient_x = numerical_partial_derivative_x(self.function_to_maximize, x, y, epsilon = epsilon)
        new_gradient_y = numerical_partial_derivative_y(self.function_to_maximize, x, y, epsilon = epsilon)
        # Perform gradient ascent
        percent_diff = 1
        i = 0
        gradient_close_to_zero = False
        while ((i < min_iterations) or (percent_diff > tolerance and i < max_iteration)):
            # Compute the partial derivatives of the function at the current point using numerical differentiation

            gradient_x = new_gradient_x
            gradient_y = new_gradient_y

            # Update the current point using the gradients and "learning rate"
            xs.append(x)
            ys.append(y)
            zs.append(z)

            newx = x + learning_rate * gradient_x
            newy = y + learning_rate * gradient_y
            newz = self.function_to_maximize(newx, newy)
            new_gradient_x = numerical_partial_derivative_x(self.function_to_maximize, newx, newy, epsilon)
            new_gradient_y = numerical_partial_derivative_y(self.function_to_maximize, newx, newy, epsilon)


            percent_diff_measurement = (abs(newz-z)/z)
            percent_diff_gradient = np.sqrt(abs(gradient_x - new_gradient_x /gradient_x)*epsilon**2 + abs(gradient_y - new_gradient_y /gradient_y)*epsilon**2)
            percent_diff = percent_diff_gradient*percent_diff_measurement*100
            x = newx
            y = newy
            z = newz

            print("Iteration {}: (x, y) = ({}, {}), f(x, y) = {}\nPercent_Diff{}".format(i + 1, x, y, self.function_to_maximize(x, y),percent_diff))
            i += 1
        if(i >= max_iteration):
            max_i = np.argmax(zs)
            x = xs[max_i]
            y = ys[max_i]
        return x, y, xs, ys, zs

    def optimize(self, range_val=180, default_steps = 3, gradient_ascent_tol=0.1, min_iterations = 10,
                 max_iterations = 20, pre_optimized = False, learning_rate = 180, gradient_epsilon = 1, warnings = True):
        """
        Optimizes the output intensity of self.measures output by moving the first two rotors instantiated in the rotator_feedback
        class.
        Uses random walk to find possible relative maxima in the specified range.
        random walks a smaller range around those points to check for more maxima. smaller range is defined as a ration
        of the initial range and the number of relative maxima that were found.
        Uses gradient ascent to reach near peak of highest relative maxima found

        :param range_val:     the value that bounds the measurement points of the rotor to (0,range_val)x(0,range,val)
                              Can be lower than 180 if maximum is within a different range
                              TODO: Add functionality for two range values if there is more symmetry on one axis
        :param default_steps: the number of different points inside the 2D-space that is bounded by the range_val that
                                the rotors will move to during the random walk
        :param gradient_ascent_tol: percentage difference between values in the gradient ascent before it terminates
                                    calculated as the 2-norm of the percentage difference for the gradients in each direction
        :param min_iterations:  minimum number of iterations before gradient_ascent terminates, even if the tolerance is
                                reached, should be a reasonably high enough value to ensure maximum finding
        :param max_iterations: Maximum number of iterations inside of gradient_ascent
        :param pre_optimized: If True, optimize skips the random walk of the total range and starts gradient ascent immediately,
                                Will find local maxima instead of global.
        :param learning_rate: a constant that is the ratio between the gradient found and the distance to next point to be measured
                                the higher the learning_rate steps the larger the step, 180 as default representing,
                                moving 90 degrees in each direction if the gradient was 0.5
                                ##TODO: Generate learning_rate inside of optimize() and gradient_ascent_2d() as function of range_value and measured values
        :param gradient_epsilon: the value in degrees that the rotors move in each direction to find the gradient
                                larger the epsilon, the lower the noise, but greater chance of finding off peak non maxima
                                values lower than 0.5 discouraged due to accuracy of rotor movement being somewhat limited
        :return: true_max_x: the position in degrees that the first rotor will be moved to
        :return: true_max_x: the value for whichever rotor is initialized first in the rotor_feedback_class
        :return: possible_z_max:the found maximimum for intensity
        :return: X, Y, Z: the set of all points that were measured during the optimization process X, and Y representing
                the positions of the first and second rotor respectivelt and Z representing the
        """
        if(range_val > 180 and warnings):
            logging.warning("range_val should be <= 180\n Set \"warnings\" parameter to True to turn off warning")
        learning_rate = learning_rate
        epsilon = gradient_epsilon
        tolerance = gradient_ascent_tol
        init_pos_1 = self._pos(rotor_num=0)
        init_pos_2 = self._pos(rotor_num=1)
        range_val = range_val
        steps = default_steps
        if not pre_optimized:
            X, Y, Z = self.move_and_measure(range=range_val, steps=steps, dry_run=False, a=init_pos_1, b=init_pos_2)

            #Estimates relative maxima from data measured so far on the order of (#of data points)^(1/2)
            order = int(np.sqrt(len(X)))
            C = argrelmax(Z, order=order)


            max_x = X[C]
            max_y = Y[C]
            max_z = Z[C]

            max_x_primes = []
            max_y_primes = []
            max_z_primes = []
            num_rel_maxes = len(C)
            for x, y, z in zip(max_x, max_y, max_z):
                print(f"Formatted pair: {x, y}")
                newX, newY, newZ = self.move_and_measure(range=range_val/num_rel_maxes, steps=steps, dry_run=False, a=x, b=y)
                X = np.append(X, newX)
                Y = np.append(Y, newY)
                Z = np.append(Z, newZ)
                rel_max = np.argmax(newZ)
                max_x_primes.append(newX[rel_max])
                max_y_primes.append(newY[rel_max])
                max_z_primes.append(newZ[rel_max])

            possible_max = np.argmax(max_z_primes)
            possible_x_max = max_x_primes[possible_max]
            possible_y_max = max_y_primes[possible_max]
        else:
            possible_x_max = init_pos_1
            possible_y_max = init_pos_2
        true_max_x, true_max_y, newXs, newYs, newZs = self.gradient_ascent_2d(learning_rate=learning_rate,
                                                                              epsilon=epsilon, tolerance=tolerance,
                                                                            max_iteration=max_iterations, min_iterations=min_iterations,
                                                                              init_x=possible_x_max, init_y=possible_y_max)
        if not pre_optimized:
            X = np.append(X, newXs)
            Y = np.append(Y, newYs)
            Z = np.append(Z, newZs)
        else:
            X, Y, Z = newXs, newYs, newZs
        self.move_to(true_max_x, 0)
        self.move_to(true_max_y, 1)
        self.wait_stop()
        possible_z_max = self.measure(measurements=1000)
        sleep(1)
        return true_max_x, true_max_y, possible_z_max, X, Y, Z

    def optimize_2(self, range_val=180, default_steps = 3, tol=0.01):
        """
        Less funnctional version of optimize(),
        finds peaks quicker but doesn't easily converge to termination
        ##TODO FIX
        """


        init_pos_1 = self._pos(rotor_num=0)
        init_pos_2 = self._pos(rotor_num=1)
        range_val = range_val
        steps = default_steps
        X, Y, Z = self.move_and_measure(range=range_val, steps=steps, dry_run=False, a=init_pos_1, b=init_pos_2)
        order = int(np.sqrt(len(X)))
        C = argrelmax(Z, order=order)
        max_x = X[C]
        max_y = Y[C]
        max_z = Z[C]

        max_x_primes = []
        max_y_primes = []
        max_z_primes = []
        num_rel_maxes = len(C)
        for x, y, z in zip(max_x, max_y, max_z):
            print(f"Formatted pair: {x, y}")
            newX, newY, newZ = self.move_and_measure(range=range_val/num_rel_maxes, steps=steps*2, dry_run=False, a=x, b=y)
            X = np.append(X, newX)
            Y = np.append(Y, newY)
            Z = np.append(Z, newZ)
            rel_max = np.argmax(newZ)
            max_x_primes.append(newX[rel_max])
            max_y_primes.append(newY[rel_max])
            max_z_primes.append(newZ[rel_max])

        possible_max = np.argmax(max_z_primes)
        possible_x_max = max_x_primes[possible_max]
        possible_y_max = max_y_primes[possible_max]
        possible_z_max = max_z_primes[possible_max]
        percent_diff = 100
        i = 1
        found_max = possible_z_max
        while(percent_diff > tol) and i <= 20:
            print(i)
            newestX, newestY, newestZ = self.move_and_measure(range=(range_val/(num_rel_maxes*(i**2))), steps=steps*2, dry_run=False, a=possible_x_max, b=possible_y_max)

            max_f = np.argmax(newestZ)
            possible_x_max = newestX[max_f]
            possible_y_max = newestY[max_f]
            found_max = newestZ[max_f]

            self.move_to(possible_x_max, 0)
            self.move_to(possible_y_max, 1)
            self.wait_stop()
            found_max = self.measure(measurements=1000)

            max_full_index = np.argmax(Z)
            full_x_max = X[max_full_index]
            full_y_max = Y[max_full_index]
            self.move_to(full_x_max, 0)
            self.move_to(full_y_max, 1)
            self.wait_stop()
            found_max_full = self.measure(measurements=1000)

            if abs((found_max_full-found_max)/found_max_full) > tol:
                i -= 1
                X = np.append(X, newestX)
                Y = np.append(Y, newestY)
                Z = np.append(Z, newestZ)
            i+=1
            percent_diff = abs((found_max_full-found_max)/found_max_full)
            print(f"Percent_Diff{percent_diff}")
        if i > 20:
            print("Fluctuations likely caused optimization to not work quickly. Check if measurement is local maximum, and re-running optimization")

            return self.optimize_2(range_val=range_val,tol = tol)
        return possible_x_max, possible_y_max, possible_z_max, X, Y, Z

    def close(self):
        self.stage[0].close()
        self.stage[1].close()
        self.daq_task.close()
        print("Closed successfully")

class RotorExperiment(EnvExperiment):

    def run(self):
        rotor1 = RotatorFeedbackChannel(ch_name="Dev1/ai0", rotator_sn=["55000741", "55105674"], dry_run=False)

        rotor1.move_to(90, 0)
        rotor1.move_to(90, 1)
        rotor1.wait_stop()
        maxX, maxY, maxZ, X, Y, Z = rotor1.optimize(range_val=180, gradient_ascent_tol=0.01, pre_optimized=False)
        fig = plt.figure()
        ax = fig.add_subplot(projection="3d")
        ax.scatter(X, Y, Z, s= 20)
        ax.scatter(maxX, maxY, maxZ, marker = "*", s = 200)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.set_title(' ')
        plt.show()
        rotor1.close()






