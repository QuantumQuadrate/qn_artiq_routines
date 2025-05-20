from time import sleep, time
from pylablib.devices import Thorlabs  # for Kinesis instrument control
import nidaqmx as daq # todo: remove
from artiq.experiment import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelmax
import logging
import os, sys

cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")
from utilities.BaseExperiment import BaseExperiment


def numerical_partial_derivative_x(function, x, y, epsilon=0.5):

    return (function(x + epsilon, y) - function(x - epsilon, y)) / (2 * epsilon)


def numerical_partial_derivative_y(function, x, y, epsilon=0.5):

    return (function(x, y + epsilon) - function(x, y - epsilon)) / (2 * epsilon)


class RotatorFeedbackChannel:

    #todo max_runs never used. delete? implement?
    def __init__(self, ch_name="Dev1/ai0", rotator_sn=['55000741','55105674'], dry_run=True,
                 max_runs=10, spc=None, rate=None, n_measurements=1, use_sampler=True, sampler_ch=0,
                 experiment=EnvExperiment):

        self.spc = spc
        self.dry_run = dry_run
        # self.daq_task = daq.Task()
        self.rate = rate

        if self.spc is None:
            self.spc = 1

        if self.rate is None or self.rate > 1e5:
            self.rate = 1e5

        self.n_measurements = n_measurements

        # self.daq_task.ai_channels.add_ai_voltage_chan(physical_channel=ch_name)

        # if self.dry_run:
        #     self.measure = None
        #     print("No dry run implemented, edit code")
        #     # use this for testing since you don't have access to artiq hardware yet

        # self.measure = self._sampler_measure
        #todo: delete
        # elif use_sampler:
        #     self.measure = self._sampler_measure
        # else:
        #     self.measure = self._nidaq_measure

        self.sampler_ch = sampler_ch
        self.experiment = experiment

        self.ser0 = rotator_sn[0]
        self.ser1 = rotator_sn[1]

        self.scl = True  # todo: whether to use physical units - this apparently has no effect

        self.stage = []
        for sn in rotator_sn:
            try:
                rotor = Thorlabs.KinesisMotor(conn=sn, scale='K10CR1', )
                self.stage.append(rotor)
                print(f"opened K10CR1 {sn}")
            except Exception as e:
                logging.error(f"Failed to connect to K10CR1 {sn}! "
                      f"\n Check that you correctly typed the SN and that it is connected.")
                print(sn, e)
                raise
        self.extra_rotors = []

        if len(rotator_sn) > 2:
            for i in range(2, len(rotator_sn)):
                try:
                    new_rotor = Thorlabs.KinesisMotor(conn=rotator_sn[i], scale='K10CR1')
                    self.extra_rotors.append(new_rotor)
                except Exception as e:
                    new_rotor.close()
                    print(sn, e)
        self.r0 = self.stage[0]
        self.r1 = self.stage[1]
        self.r0.set_position_reference()
        self.r1.set_position_reference()

        #(f"{self._pos(0)} : {self._pos(1)}; intensity:{np.sum(self.daq_task.read(number_of_samples_per_channel=1))}")

    # todo delete unused args
    def print_connections(self, all = True, rotors = False, task = False):
        if rotors is True:
            for r in self.stage:
                print(r)
            for r in self.extra_rotors:
                print(r)

    def print_pos(self, rotor_num=0):
        print("position = ", self._pos(rotor_num))

    def wait_stop(self, rotor_num=2):
        if rotor_num == 2:
            self.stage[0].wait_for_stop()
            self.stage[1].wait_for_stop()
        else:
            self.stage[rotor_num].wait_for_stop()

    def is_moving(self, rotor_num=None):
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

    # todo: delete. nidaq resolution is too low
    # def _nidaq_measure(self, measurements=None):
    #     """Read photodetector with USB-NiDaq"""
    #     if measurements is None:
    #         measurements = self.spc
    #
    #     data = self.daq_task.read(number_of_samples_per_channel=measurements)
    #
    #     return np.mean(data)

    @kernel
    def measure(self, measurements: TInt32 = 10, averaging_delay: TFloat = 1*ms) -> TFloat:
        """Read photodetector with Sinara Sampler card 1"""

        buffer = [0.0]*8
        data = [0.0]*measurements
        for i in range(measurements):
            delay(averaging_delay)
            self.experiment.sampler1.sample(buffer)
            data[i] = buffer[self.sampler_ch]

        return np.mean(data)

    def q_h_gen(self, motor_range, steps, center_angle1, center_angle2, random=False):
        """

        :param motor_range: width and height of the range that the random/pre-generated points will be within
        :param steps: amount of points
        :param center_angle1: center of axis of the range for "rotor1"
        :param center_angle2: center of axis of the range for "rotor2"
        :param random: if False generates purely uniformly random points with specified motor_range and center
        :return: x, y arrays of the random points where len(x) = len(y) = steps
        """
        if not random:
            motor_range /= 2
            rads = np.arange(0, 2 * np.pi, (2 * np.pi) / steps)
            x = motor_range * (np.cos(rads)) + center_angle1
            y = motor_range * (np.sin(rads)) + center_angle2
            x1 = motor_range / 2 * (np.cos(rads)) + center_angle1
            y1 = motor_range / 2 * (np.sin(rads)) + center_angle2
            x2 = motor_range / 4 * (np.cos(rads)) + center_angle1
            y2 = motor_range / 4 * (np.sin(rads)) + center_angle2
            x0 = np.append(x, x1)
            x0 = np.append(x0, x2)
            y0 = np.append(y, y1)
            y0 = np.append(y0, y2)

        if random is True:
            x0 = motor_range*(np.random.rand(steps) - 1/2) + center_angle1
            y0 = motor_range*(np.random.rand(steps) - 1/2) + center_angle2

        return x0, y0

    def move_and_measure(self, motor_range=None, steps=None, rotor1_initial_degrees=0, rotor2_initial_degrees=0,
                         dry_run=True, warnings=True):
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
        :param: motor_range: motor_range of values in both directions that will be walked through
        :param: steps: the number of 2d points that will be moved between during the random walk
        :param: a = initial position of rotor1
        :param: b = initial position of rotor2
        :return: returns arrays waveplate1_angles_deg, waveplate2_angles_deg, and measurements,
            which are the arrays of waveplate angles in degrees and the corresponding measured voltages
        """
        ang1, ang2 = self.q_h_gen(steps=steps, motor_range=motor_range,
                                  center_angle1=rotor1_initial_degrees,
                                  center_angle2=rotor2_initial_degrees, random=True)
        m_ang1 = []
        m_ang2 = []
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

                # measures while rotors are in motion
                # while self.is_moving():
                #     m_ang1 = np.append(m_ang1, self._pos(rotor_num=0))
                #     m_ang2 = np.append(m_ang2, self._pos(rotor_num=1))
                #     measurement = self.measure(measurements=self.n_measurements)
                #     measurements = np.append(measurements, measurement)
                #
                #     i += 1
                #     sleep(0.3)
                # self.stage[1].wait_for_stop()
                # self.stage[0].wait_for_stop()

                self.stage[1].wait_for_stop()
                self.stage[0].wait_for_stop()
                measurement = self.measure(measurements=self.n_measurements)
                measurements = np.append(measurements, measurement)

        waveplate1_angles_deg = m_ang1 % 180
        waveplate2_angles_deg = m_ang2 % 180

        return waveplate1_angles_deg, waveplate2_angles_deg, measurements

    def move_by(self, degrees, rotor_num=-1, velocity=None, r1=None):
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

    def move_to(self, degrees, rotor_num, velocity=None, r1=None):
        r = r1
        if rotor_num == -1 and r is None:
            logging.warning("Must provide rotor number or instantiated rotor class: no movement will be made")
            return -1
        else:
            r = self.stage[rotor_num]

        r.move_to(degrees)
        r.wait_move()
        return 1

    def function_to_maximize(self, x, y):
        """

        :param x: the position that the first initialized rotor "self.stage[0]" will be moved to
        :param y: the position that the first initialized rotor "self.stage[0]" will be moved to
        :return: measurement value read for NiDaq
        """
        self.move_to(x, 0)
        self.move_to(y, 1)

        # ensures that the rotors are at the correct position before any measurement is made
        self.wait_stop(rotor_num=0)
        self.wait_stop(rotor_num=1)
        measurement = self.measure(measurements=1000)
        sleep(0.1)  # sleeping to ensure that the measurement is completed before the rotors start moving again
        return measurement

    # todo: use max runs arg from class instance instead of kwarg?
    def gradient_ascent_2d(self, learning_rate, tolerance, min_iterations=10, max_iteration=100, init_angle1=0, init_angle2=0,
                           epsilon=0.5):
        # Start with initial guesses for the maximum (x, y) # are x/y the angles?
        x = init_angle1
        y = init_angle2
        z = self.function_to_maximize(x, y)
        xs = []
        ys = []
        zs = []
        new_gradient_angle1 = numerical_partial_derivative_x(self.function_to_maximize, x, y, epsilon=epsilon)
        new_gradient_angle2 = numerical_partial_derivative_y(self.function_to_maximize, x, y, epsilon=epsilon)
        # Perform gradient ascent
        percent_diff = 1
        i = 0

        # todo: unused, delete or implement
        gradient_close_to_zero = False
        while (i < min_iterations) or (percent_diff > tolerance and i < max_iteration):
            # Compute the partial derivatives of the function at the current point using numerical differentiation

            gradient_angle1 = new_gradient_angle1
            gradient_angle2 = new_gradient_angle2

            # Update the current point using the gradients and "learning rate"
            xs.append(x)
            ys.append(y)
            zs.append(z)

            newx = x + learning_rate * gradient_angle1
            newy = y + learning_rate * gradient_angle2
            newz = self.function_to_maximize(newx, newy)
            new_gradient_angle1 = numerical_partial_derivative_x(self.function_to_maximize, newx, newy, epsilon)
            new_gradient_angle2 = numerical_partial_derivative_y(self.function_to_maximize, newx, newy, epsilon)

            percent_diff_measurement = (abs(newz-z)/z)
            percent_diff_gradient = np.sqrt(
                abs(gradient_angle1 - new_gradient_angle1 /gradient_angle1)*epsilon**2
                + abs(gradient_angle2 - new_gradient_angle2 /gradient_angle2)*epsilon**2)
            percent_diff = percent_diff_gradient*percent_diff_measurement*100
            x = newx
            y = newy
            z = newz

            print("Iteration {}: (x, y) = ({}, {}), f(x, y) = {}\nPercent_Diff{}"
                  .format(i + 1, x, y, self.function_to_maximize(x, y),percent_diff))
            i += 1
        if i >= max_iteration:
            max_i = np.argmax(zs)
            x = xs[max_i]
            y = ys[max_i]
        return x, y, xs, ys, zs

    def optimize(self, range_val=180, default_steps=3, gradient_ascent_tol=0.1, min_iterations=10,
                 max_iterations=20, pre_optimized=False, learning_rate=180, gradient_epsilon=1, warnings=True):
        """
        Optimizes the output intensity of self.measure output by moving the first two rotors instantiated in the rotator_feedback
        class.
        Uses random walk to find possible relative maxima in the specified range.
        random walks a smaller range around those points to check for more maxima. smaller range is defined as a ratio
        of the initial range and the number of relative maxima that were found.
        Uses gradient ascent to reach near peak of highest relative maxima found

        :param range_val:     the value that bounds the measurement points of the rotor to (0,range_val)angle1(0,range,val)
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
        :return: true_max_angle1: the position in degrees that the first rotor will be moved to
        :return: true_max_angle2: the value for whichever rotor is initialized first in the rotor_feedback_class
        :return: possible_z_max:the found maximimum for intensity
        :return: waveplate1_angles, waveplate2_angles, Z: the set of all points that were measured during the optimization process waveplate1_angles, and waveplate2_angles representing
                the positions of the first and second rotor respectively and Z representing the
        """
        if range_val > 180 and warnings:
            logging.warning("range_val should be <= 180\n Set \"warnings\" parameter to False to turn off warning")
        epsilon = gradient_epsilon
        tolerance = gradient_ascent_tol
        init_pos_1 = self._pos(rotor_num=0)
        init_pos_2 = self._pos(rotor_num=1)
        steps = default_steps
        if not pre_optimized:
            waveplate1_angles, waveplate2_angles, measurements = self.move_and_measure(motor_range=range_val,
                                                                                       steps=steps,
                                                                                       rotor1_initial_degrees=init_pos_1,
                                                                                       rotor2_initial_degrees=init_pos_2,
                                                                                       dry_run=False)

            # Estimates relative maxima from data measured on the order of (#of data points)^(1/2)
            order = int(np.sqrt(len(waveplate1_angles)))
            relative_maxima_indices = argrelmax(measurements, order=order) # todo more useful variable name

            relative_maxima_angles1 = waveplate1_angles[relative_maxima_indices]
            relative_maxima_angles2 = waveplate2_angles[relative_maxima_indices]
            relative_maxima_measurements = measurements[relative_maxima_indices]

            # todo: what is meant by primes?
            max_angle1_primes = []
            max_angle2_primes = []
            max_z_primes = []
            num_rel_maxes = len(relative_maxima_indices)

            # this loop does a kind of 1-layer deep tree search using the previously found relative maxima
            for angle1, angle2, measured_local_max in zip(relative_maxima_angles1, relative_maxima_angles2,
                                                          relative_maxima_measurements):
                print(f"Formatted pair: {angle1, angle2}")
                new_angles1, new_angles2, new_measurements = self.move_and_measure(
                    motor_range=range_val / num_rel_maxes, steps=steps, rotor1_initial_degrees=angle1,
                    rotor2_initial_degrees=angle2, dry_run=False)
                waveplate1_angles = np.append(waveplate1_angles, new_angles1)
                waveplate2_angles = np.append(waveplate2_angles, new_angles2)
                measurements = np.append(measurements, new_measurements)
                rel_max = np.argmax(new_measurements)
                max_angle1_primes.append(new_angles1[rel_max])
                max_angle2_primes.append(new_angles2[rel_max])
                max_z_primes.append(new_measurements[rel_max])

            possible_max = np.argmax(max_z_primes)
            possible_angle1_max = max_angle1_primes[possible_max]
            possible_angle2_max = max_angle2_primes[possible_max]
        else:
            possible_angle1_max = init_pos_1
            possible_angle2_max = init_pos_2
        true_max_angle1, true_max_angle2, newXs, newYs, newZs = self.gradient_ascent_2d(
            learning_rate=learning_rate,
            epsilon=epsilon, tolerance=tolerance,
            max_iteration=max_iterations, min_iterations=min_iterations,
            init_angle1=possible_angle1_max, init_angle2=possible_angle2_max)
        if not pre_optimized:
            waveplate1_angles = np.append(waveplate1_angles, newXs)
            waveplate2_angles = np.append(waveplate2_angles, newYs)
            measurements = np.append(measurements, newZs)
        else:
            waveplate1_angles, waveplate2_angles, measurements = newXs, newYs, newZs
        self.move_to(true_max_angle1, 0)
        self.move_to(true_max_angle2, 1)
        self.wait_stop()
        possible_z_max = self.measure(measurements=1000)
        sleep(1)
        return true_max_angle1, true_max_angle2, possible_z_max, waveplate1_angles, waveplate2_angles, measurements

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
        X, Y, Z = self.move_and_measure(motor_range=range_val, steps=steps, rotor1_initial_degrees=init_pos_1,
                                        rotor2_initial_degrees=init_pos_2, dry_run=False)
        order = int(np.sqrt(len(X)))
        C = argrelmax(Z, order=order)
        max_angle1 = X[C]
        max_angle2 = Y[C]
        max_z = Z[C]

        max_angle1_primes = []
        max_angle2_primes = []
        max_z_primes = []
        num_rel_maxes = len(C)
        for x, y, z in zip(max_angle1, max_angle2, max_z):
            print(f"Formatted pair: {x, y}")
            newX, newY, newZ = self.move_and_measure(motor_range=range_val / num_rel_maxes, steps=steps * 2,
                                                     rotor1_initial_degrees=x, rotor2_initial_degrees=y, dry_run=False)
            X = np.append(X, newX)
            Y = np.append(Y, newY)
            Z = np.append(Z, newZ)
            rel_max = np.argmax(newZ)
            max_angle1_primes.append(newX[rel_max])
            max_angle2_primes.append(newY[rel_max])
            max_z_primes.append(newZ[rel_max])

        possible_max = np.argmax(max_z_primes)
        possible_angle1_max = max_angle1_primes[possible_max]
        possible_angle2_max = max_angle2_primes[possible_max]
        possible_z_max = max_z_primes[possible_max]
        percent_diff = 100
        i = 1
        found_max = possible_z_max
        while(percent_diff > tol) and i <= 20:
            print(i)
            newestX, newestY, newestZ = self.move_and_measure(motor_range=(range_val / (num_rel_maxes * (i ** 2))),
                                                              steps=steps * 2, rotor1_initial_degrees=possible_angle1_max,
                                                              rotor2_initial_degrees=possible_angle2_max, dry_run=False)

            max_f = np.argmax(newestZ)
            possible_angle1_max = newestX[max_f]
            possible_angle2_max = newestY[max_f]
            found_max = newestZ[max_f]

            self.move_to(possible_angle1_max, 0)
            self.move_to(possible_angle2_max, 1)
            self.wait_stop()
            found_max = self.measure(measurements=1000)

            max_full_index = np.argmax(Z)
            full_angle1_max = X[max_full_index]
            full_angle2_max = Y[max_full_index]
            self.move_to(full_angle1_max, 0)
            self.move_to(full_angle2_max, 1)
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
        return possible_angle1_max, possible_angle2_max, possible_z_max, X, Y, Z

    # def fastimize(self):
    #     """
    #     For optimizing when we are close to optimium polarization, which is typically the case
    #     :return:
    #     """



    def close(self):
        self.stage[0].close()
        self.stage[1].close()
        # self.daq_task.close()
        print("Closed successfully")


# class RotorExperiment(EnvExperiment):
#
#     def run(self):
#         rotor1 = RotatorFeedbackChannel(ch_name="Dev1/ai0", rotator_sn=["55000741", "55105674"], dry_run=False)
#
#         rotor1.move_to(90, 0)
#         rotor1.move_to(90, 1)
#         rotor1.wait_stop()
#         maxX, maxY, maxZ, X, Y, Z = rotor1.optimize(range_val=180, gradient_ascent_tol=0.01, pre_optimized=False)
#         fig = plt.figure()
#         ax = fig.add_subplot(projection="3d")
#         ax.scatter(X, Y, Z, s=20)
#         ax.scatter(maxX, maxY, maxZ, marker="*", s=200)
#         ax.set_angle1label('x')
#         ax.set_ylabel('y')
#         ax.set_zlabel('z')
#         ax.set_title(' ')
#         plt.show()
#         rotor1.close()


class FORTPolarizationOptimizer(EnvExperiment):

    def build(self):
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        # does it matter which waveplate is which?
        self.setattr_argument('K10CR1_serial_numbers', StringValue('["55000759", "55000740"]'))
        self.setattr_argument('optimize_kwargs',
                              StringValue("{'range_val':180, 'gradient_ascent_tol':0.01, 'pre_optimized':False}"))

        self.base.set_datasets_from_gui_args()

    def prepare(self):
        self.base.prepare()
        self.k10cr1_SNs = eval(self.K10CR1_serial_numbers)
        assert len(self.k10cr1_SNs) == 2, "there should be two serial numbers"
        self.laser_stabilizer.dds_names = ['dds_FORT']  # the only AOM to which we're feeding back

        assert hasattr(self, 'core'), "uh oh"
        self.rotors = RotatorFeedbackChannel(rotator_sn=self.k10cr1_SNs,
                                             dry_run=False, sampler_ch=0, experiment=self)

    @kernel
    def sinara_hardware_setup(self):
        """
        Initialization and turning stuff on.
        :return:
        """
        self.base.initialize_hardware()

        self.core.reset()

        self.dds_FORT.sw.on()
        delay(100*ms)
        self.laser_stabilizer.run()
        self.dds_FORT.sw.on()

        logging.info("Sinara hardware setup - done")

    def run(self):

        self.sinara_hardware_setup()

        start = time()
        result = self.rotors.optimize(**eval(self.optimize_kwargs))
        max_angle1, max_angle2, max_measurement, angles1, angles2, measurements = result
        logging.info(f"rotor optimization completed in {time()-start:.2f}s")
        fig = plt.figure()
        ax = fig.add_subplot(projection="3d")
        ax.scatter(angles1, angles2, measurements, s=20)
        ax.scatter(max_angle1, max_angle2, max_measurement, marker="*", s=200)
        ax.set_xlabel('angle 1 (deg.)')
        ax.set_ylabel('angle 2 (deg.)')
        ax.set_zlabel('Photodetector signal (arb.)')
        ax.set_title('FORT optimization:'+str(self.scheduler.rid))
        plt.show()
        self.rotors.close()


if __name__ == '__main__':
    pass






