"""
WORK IN PROGRESS: uses measurement after a fiber and PBS to estimate stokes parameters and the polarizing action of a optical fiber as a
result
"""
from time import sleep, time
from pylablib.devices import Thorlabs # for Kinesis instrument control
import nidaqmx as daq
from artiq.experiment import *
import nidaqmx.constants as daq_constants
from nidaqmx.errors import DaqError, DaqWarning
from nidaqmx.error_codes import DAQmxErrors, DAQmxWarnings
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from ArbitraryRetarder import *
from scipy.optimize import minimize
from scipy.stats import norm
# a dictionary specifying channels we want to optimize

class RotatorFeedbackChannel():

    def __init__(self,ch_name = "Dev1/ai0", dds_channel=0, rotator_sn=['55000741','55105674'], dry_run = True ,
                 max_runs=10, leave_laser_on=False, spc = None, rate = None, plate_config = None):

        self.spc = spc
        self.dry_run = dry_run
        self.daq_task = daq.Task()
        self.rate = rate

        if self.spc is None:
            self.spc = 1

        if self.rate is None or self.rate > 1e5:
            self.rate = 1e5

        self.measure_sim = plate_config_measure(plate_config)

        self.daq_task.ai_channels.add_ai_voltage_chan(physical_channel=ch_name)

        if self.dry_run:
            self.measure = self.measure_sim  # use this for testing since you don't have access to artiq hardware yet
        else:
            self.measure = self._measure

        self.ser0 = rotator_sn[0]
        self.ser1 = rotator_sn[1]

        self.scl = True  # whether to use physical units - this apparently has no effect
        self.stage = [Thorlabs.KinesisMotor(conn=self.ser0, scale='K10CR1', ),
                      Thorlabs.KinesisMotor(conn=self.ser1, scale='K10CR1', )]

        self.r0 = self.stage[0]
        self.r1 = self.stage[1]
        self.r0.set_position_reference()
        self.r1.set_position_reference()

        print(f"{self._pos(0)} : {self._pos(1)}; intensity:{np.sum(self.daq_task.read(number_of_samples_per_channel=1))}")

    def expected_improvement(self, x, ang1, ang2, data, xi=0.01):
        arb_angs = (x[0])
        E = x[3]
        background = x[4]
        ang1-=x[5]
        ang2-=x[6]
        angs = [ang1, ang2, arb_angs]
        mean_intensity = self.measure_sim(angles=angs, E = E, background=background)  # Replace with your actual prediction method
        std_intensity = np.sqrt(mean_intensity)  # Assuming you have a method for predicting standard deviation
        if data is not None:
            best_observed = np.max(data[2])# Assuming the best observed intensity is the maximum
        else:
            best_observed = np.random()*5+2
        z = (mean_intensity - best_observed - xi) / std_intensity
        ei = (mean_intensity - best_observed - xi) * norm.cdf(z) + std_intensity * norm.pdf(z)
        return ei

    def select_next_point(self, data, x0, xi = 0.01, points = 5):
        ang1new = np.array([])
        ang2new = np.array([])
        for i in range(0, points):
            ang1_random = np.random.uniform(0,180, 100)
            ang2_random = np.random.uniform(0,180, 100)
            ei = np.array([self.expected_improvement(x0, a1, a2, data) for a1, a2 in zip(ang1_random, ang2_random)])
            ang1new = np.append(ang1new, ang1_random[np.argmax(ei)])
            ang2new = np.append(ang2new, ang2_random[np.argmax(ei)])
        return ang1new, ang2new

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


    def q_h_gen(self, xy_range, steps, center_x, center_y, random=False):
        """
        Return a list of x coordinates and y coordinates

        What are x and y in the lab?

        Parameters:
        ----------
        xy_range : float
            The maximum distance from the center for the generated points.
        steps : int
            The number of points to generate.
        center_x : float
            The x-coordinate of the center point.
        center_y : float
            The y-coordinate of the center point.
        random : bool, optional
            If True, generate points randomly within the range. Default is False.

        Returns:
        -------
        tuple of numpy.ndarray
            Two arrays containing the x and y coordinates of the generated points.
        """
        if not random:
            xy_range /= 2
            rads = np.arange(0, 2 * np.pi, (2 * np.pi) / steps)
            x = xy_range * (np.cos(rads)) + center_x
            y = xy_range * (np.sin(rads)) + center_y
            x1 = xy_range / 2 * (np.cos(rads)) + center_x
            y1 = xy_range / 2 * (np.sin(rads)) + center_y
            x2 = xy_range / 4 * (np.cos(rads)) + center_x
            y2 = xy_range / 4 * (np.sin(rads)) + center_y
            x0 = np.append(x, x1)
            x0 = np.append(x0, x2)
            y0 = np.append(y, y1)
            y0 = np.append(y0, y2)
        if random is True:
            x0 = np.random.rand(steps) * xy_range + center_x
            y0 = np.random.rand(steps) * xy_range + center_y

        return x0, y0

    def move_and_measure(self, del_phi=None, angular_range=None, steps=None, a=0, b=0, theta_1=0, theta_2=0,
                         dry_run=True, E=None, start_time=None, points=None):
        """
        Using pre selected point/randomly generated ones, the rotors are moved between points and then stopped.
        While the rotors are in the process of moving, we are keeping track of the angle it is at,
        and then the measurement. Due to time to measure, acceleration and speed of rotor, there is slight discrepancy
        between measured angle and actual angle for some given intensity measurement, but the non acceleration component
        of the discrepancy is solved for in the optimize portion of the code, and the acceleration error component
        seems to be negligible, when the goal is achieving maximum intensity.
        """
        if points is None:
            ang1, ang2 = self.q_h_gen(xy_range=angular_range, steps=steps, center_x=a, center_y=b, random=True)
            ang1 -= theta_1
            ang2 -= theta_2
        else:
            ang1, ang2 = points
        m_ang1 = []
        m_ang2 = []
        measure = self.measure
        measurements = []
        if dry_run or self is None:
            for h, q in zip(ang1, ang2):
                # todo delete
                # measurements.append(np.sum(measure(q_ang=q, h_ang=h, del_phi=del_phi, E=E))) # measure doesn't take any of these arguments...
                measurements.append(np.sum(measure()))
        else:
            for a1, a2 in zip(ang1, ang2):
                self.r0.move_to(a1)
                self.r1.move_to(a2)
                i = 0
                sleep(1)
                while self.is_moving():
                    m_ang1 = np.append(m_ang1, self._pos(rotor_num=0))
                    m_ang2 = np.append(m_ang2, self._pos(rotor_num=1))
                    measurement = measure(measurements=1)
                    measurements = np.append(measurements, measurement)
                    # print("intermediate measurement")
                    # print(q_ang)
                    i += 1
                    sleep(0.05)
                self.stage[1].wait_for_stop()
                self.stage[0].wait_for_stop()

        return m_ang1%180, m_ang2%180, measurements

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

    @staticmethod
    def max_with_tolerance(self, h, q, m, peak, tol):
        m_indices = [(peak * (1 + tol)) > i > (peak * (1 - tol)) for i in m] # this is never used... # todo delete
        for i, j, m_1 in zip(h, q, m):
            if m_1 == max(m) and ((peak * (1 + tol)) > m_1 > (peak * 1 - tol)):
                return i, j, m_1
        return None, None, None

    def move_to(self, degrees, rotor_num = 0, velocity = None, r1=None):
        r = r1
        if rotor_num == -1 and r is None:
            r = self.stage[0]
        else:
            r = self.stage[rotor_num]

        r.move_to(degrees)
        r.wait_move()
        return 1

    def objective_func(self, x, ang1_data, ang2_data, mdata):
        del_phi = x[0]
        E = x[1]
        a = x[2]
        p1 = x[3]
        p2 = x[4]
        measurements = mdata
        error = 0
        arb_r_angs = [del_phi]
        for ang1, ang2, m in zip(ang1_data, ang2_data, measurements):
            angs = [ang1-p1, ang2-p2, arb_r_angs]
            error += (np.sum(np.array(self.measure_sim(angs, E = E, background = a))-m)) ** 2
        return error

    def optimize(self, data=None, x0=None, bounds=None, cons=None, range_val=180,
                 terminate=False, tol=0.2, steps = 5):
        """
        :param data:
        :param x0:
        :param bounds:
        :param cons:
        :param range_val:
        :param terminate:
        :param alpha:
        :param beta:
        :param tol:
        :return:
        """
        def constraint(x, args):
            """

            :param x[1]: guessed parameter value for amplitude
            :param x[2]: guessed paramter value for background amplitude
            :param args: maximum measured value from data
            return value multiplied by 1e10 so that the scip.py minimize weights the peak vbeing correct over the indivi
            -dual point matching.
            :return: E+a = peak;  (peak-max_val) > 0
            """
            return (x[1] + x[2] - args)

        cons = ({'type': 'ineq',
                 'fun': constraint,
                 'args': (max(data[2]),)
                 },)
        r0 =  self.r0
        r1 =  self.r1
        method = 'trust-constr'
        trust_constr_opts = {'disp': True, 'barrier_tol': 1e-8, 'xtol': 1e-15,
                             'initial_constr_penalty': 1, 'maxiter': 5e3, }
        if data is not None:
            pot_max = max(data[2])
            pot_min = min(data[2])
            bounds = [(0, pi), (pot_max - pot_min, np.inf),
                      (0, pot_max), (0, 180), (0,180)]
        else:
            pot_max = rand() * 3
            pot_min = 0
        if x0 is None:
            x0 = [rand() * pi, rand() * pi, rand() * pi, pot_max - pot_min, 0, 0, 0]
            sleep(1)
        result = minimize(self.objective_func, x0=x0, args=data, bounds=bounds, method=method,
                          constraints=cons, options=trust_constr_opts)
        print(result.x)
        x = result.x
        del_phi, E, a, p1, p2 = x

        ang1 = np.linspace(0, range_val, 30)
        ang2 = np.linspace(0, range_val, 30)

        X, Y = np.meshgrid(ang1, ang2)

        arb_angs = del_phi
        angs = [X, Y, arb_angs]

        Z1 = self.measure_sim(angles=angs, E=E, background=a)
        peak = (E + a)

        c = np.argmax(Z1)

        maxX, maxY, maxZ = X[c // 30][c % 30], Y[c // 30][c % 30], Z1[c // 30][c % 30]

        self.r0.move_to(maxX)
        self.r1.move_to(maxY)
        self.r0.wait_for_stop()
        self.r1.wait_for_stop()
        test_measurement = self.measure(measurements=1000)
        sleep(1)
        percent_diff = (abs(test_measurement - peak) / peak)
        within_range = percent_diff < tol
        terminate = within_range or terminate
        print(percent_diff)
        if not terminate:
            x = del_phi, E, a, p1, p2
            gen_points = self.select_next_point(data, x0=x, points=steps)
            newh, newq, newm = self.move_and_measure(angular_range=180, steps=5, a=0, b=0, dry_run=False,
                                                     points=gen_points)
            h, q, m = data
            h = np.append(h, newh)
            q = np.append(q, newq)
            m = np.append(m, newm)
            data = h, q, m
            # Before trying to reoptimize using new data it checks if any of the new points fall within the tolerance of
            # previous peak found by scipy.minimize
            checkMax1, checkMax2, checkMaxM = self.max_with_tolerance(*data, peak, tol)
            if checkMax1 is not None and checkMax1 is not None:
                r0.move_to(checkMax1)
                r1.move_to(checkMax2)
                r0.wait_for_stop()
                r1.wait_for_stop()
                test_measurement = self.measure(measurements=1000)
                # Double checks to ensure the point measured to have a max within the tolerance has that max again
                within_range = test_measurement >= peak * (1 - tol) and test_measurement <= peak * (1 + tol)
                print(within_range)
                print(abs(test_measurement-peak)/(peak * (1 - tol)))
                if within_range:
                    print("returning")
                    return X, Y, Z1, checkMax1, checkMax2, checkMaxM, data
            # if no peak is found within the total new range, it is re optimized to find a new peak to search for
            return self.optimize(x0=x0, data=data, terminate=False, tol=tol)
        else:
            return X, Y, Z1, maxX, maxY, maxZ, data

    def close(self):
        self.stage[0].close()
        self.stage[1].close()
        self.daq_task.close()
        print("Closed successfully")

class RotorExperiment2(EnvExperiment):

    def run(self):
        rotor1 = RotatorFeedbackChannel(ch_name="Dev1/ai0", rotator_sn=["55000741", "55105674"], dry_run=False, plate_config = (qwp, hwp, fiber))
        init_pos_1 = rotor1._pos(rotor_num=0)
        init_pos_2 = rotor1._pos(rotor_num=1)
        range_val = 180
        steps = 15
        del_phi = None
        E_0 = None
        data = rotor1.move_and_measure(del_phi=del_phi, angular_range=range_val, steps=steps, a=init_pos_1,
                                       b=init_pos_2, dry_run=False, E=E_0)

        X, Y, Z1, maxX, maxY, maxZ, data = rotor1.optimize(data=data, tol=0.03)
        fig = plt.figure()

        ax = fig.add_subplot(projection="3d")

        ax.scatter(maxX, maxY, maxZ, marker="*", s=20)
        ax.scatter(X, Y, Z1, marker=".", s=5)

        ax.scatter(data[0], data[1], data[2])

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.set_title(' ')
        plt.show()
        rotor1.close()






