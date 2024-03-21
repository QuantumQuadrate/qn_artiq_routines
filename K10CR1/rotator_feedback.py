
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
# a dictionary specifying channels we want to optimize

class RotatorFeedbackChannel():

    def __init__(self,ch_name = "Dev1/ai0", dds_channel=0, rotator_sn=['55105674', '55000741'], dry_run = True ,
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


    """
    Self explanatory, generated points for two wave plates, named q and h here, represented by y0 and x0
    """

    def q_h_gen(self, range, steps, center_x, center_y, random=False):

        if not random:
            range /= 2
            rads = np.arange(0, 2 * np.pi, (2 * np.pi) / steps)
            x = range * (np.cos(rads)) - center_x
            y = range * (np.sin(rads)) - center_y
            x1 = range / 2 * (np.cos(rads)) - center_x
            y1 = range / 2 * (np.sin(rads)) - center_y
            x2 = range / 4 * (np.cos(rads)) - center_x
            y2 = range / 4 * (np.sin(rads)) - center_y
            x0 = np.append(x, x1)
            x0 = np.append(x0, x2)
            y0 = np.append(y, y1)
            y0 = np.append(y0, y2)
        if random is True:
            x0, x1, x2 = np.random.rand(steps) * range, np.random.rand(steps) * range / 2, np.random.rand(steps) * range / 4
            x0 = np.append(x0, x1)
            x0 = np.append(x0, x2)
            y0, y1, y2 = np.random.rand(steps) * range, np.random.rand(steps) * range / 2, np.random.rand(steps) * range / 4
            y0 = np.append(y0, y1)
            y0 = np.append(y0, y2)

        return x0, y0
    """
    Using pre selected point/randomly generated ones, the rotors are moved between points and then stopped.
    While the rotors are in the process of moving, we are keeping track of the angle it is at, and then the measurement.
    Due to time to measure, acceleration and speed of rotor, there is slight discrepancy between measured angle and 
    actual angle for some given intensity measurement, but the non acceleration component of the discrepancy is solved 
    for in the optimize portion of the code, and the acceleration error component seems to be negligible, when the goal
    is achieving maximum intensity.
    """
    def move_and_measure(self, theta=None, eta=None, phi=None, range=None, steps=None, a=0, b=0, theta_h=0, theta_q=0,
                         dry_run=True, E=None, background = 0):
        h_ang, q_ang = self.q_h_gen(steps=steps, range=range, center_x=-a, center_y=-b, random=True)
        q_ang -= theta_q
        h_ang -= theta_h
        mq_ang = []
        mh_ang = []
        # q_ang.sort()
        # h_ang.sort()
        measure = self.measure
        measurements = []
        if dry_run or self is None:
            arb_angles = [theta, eta, phi]
            for q, h in zip(q_ang, h_ang):
                angles = [q, h, arb_angles]
                measurements.append(np.sum(self.measure(angles, E = E, background = background)))
        else:
            for q, h in zip(q_ang, h_ang):
                self.r1.move_to(q)
                self.r0.move_to(h)
                i = 0
                sleep(0.3)
                while self.is_moving():
                    mq_ang = np.append(mq_ang, self._pos(rotor_num=1))
                    mh_ang = np.append(mh_ang, self._pos(rotor_num=0))
                    measurement = measure(measurements=1)
                    measurements = np.append(measurements, measurement)
                    # print("intermediate measurement")
                    # print(q_ang)
                    i += 1
                    if not i % 10:
                        print("still reading...")
                    sleep(0.3)
                self.stage[1].wait_for_stop()
                self.stage[0].wait_for_stop()

        return mh_ang, mq_ang, measurements
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

    def max_with_tolerance(self, h, q, m, peak, tol):
        m_indices = [i < (peak * (1 + tol)) and i > (peak * (1 - tol)) for i in m]
        m_indices = [j == max(m) and t for j, t in zip(m, m_indices)]
        for i, j, m in zip(h, q, m):
            if i != -361 and j != -361:
                return i, j, m
        return None, None, None

    def move_to(self, degrees, rotor_num = 0, velocity = None, r1 = None):
        r = r1
        if rotor_num == -1 and r is None:
            r = self.stage[0]
        else:
            r = self.stage[rotor_num]

        r.move_to(degrees)
        r.wait_move()
        return 1

    def objective_func(self, x, ang1_data, ang2_data, mdata):
        theta = x[0]
        eta = x[1]
        phi = x[2]
        E = x[3]
        phase_1 = x[4]
        phase_2 = x[5]
        a = x[6]
        measurements = mdata
        error = 0
        arb_r_angs = (theta,eta,phi)
        for ang1, ang2, m in zip(ang1_data, ang2_data, measurements):
            angs = (ang1-phase_1, ang2-phase_2, arb_r_angs)
            error += (np.sum(np.array(self.measure_sim(angs, E = E, background = a))-m)) ** 2

    def optimize(self,
                 data=None, x0=None, bounds=None, cons=None, range_val=180,
                 terminate=False, alpha=0, beta=0, tol=0.2):
        cons = ({'type': 'ineq',
                 'fun': constraint,
                 'args': (max(data[2]),)
                 },)
        r0 =  self.r0
        r1 =  self.r1
        method = 'trust-constr'
        trust_constr_opts = {'disp': True, 'barrier_tol': 1e-8, 'xtol': 1e-12,
                             'initial_constr_penalty': 1, 'maxiter': 5e3, }
        if data is not None:
            pot_max = max(data[2])
            pot_min = min(data[2])
            bounds = [(0, pi), (0, pi), (0, pi), (pot_max - pot_min, np.inf), (0, pi), (0, pi),
                      (-0.5, pot_max)]
        else:
            pot_max = rand() * 3
        if x0 is None:
            x0 = [rand() * pi, rand() * pi, rand() * pi, pot_max - pot_min, 0, 0, 0]
            sleep(1)
        result = minimize(self.objective_func, x0=x0, args=data, bounds=bounds, method=method,
                          constraints=cons, options=trust_constr_opts)

        x = result.x
        theta, eta, phi, E, p_q, p_h, a = x

        h = np.linspace(-range_val / 2 + alpha, range_val / 2 + alpha, 30)
        q = np.linspace(-range_val / 2 + beta, range_val / 2 + beta, 30)

        X, Y = np.meshgrid(h, q)

        Z1 = measure(q_ang=Y, h_ang=X, theta=theta, eta=eta, phi=phi, theta_q=p_q, theta_h=p_h, E=E, a=a)
        peak = (E + a)

        c = np.argmax(Z1)

        maxX, maxY, maxZ = X[c // 30][c % 30], Y[c // 30][c % 30], Z1[c // 30][c % 30]

        self.r0.move_to(maxX)
        self.r1.move_to(maxY)
        self.r0.wait_for_stop()
        self.r1.wait_for_stop()
        test_measurement = self.measure(measurements=1000)
        sleep(1)
        print((abs(test_measurement - peak) / peak))
        within_range = test_measurement >= peak * (1 - tol) and test_measurement <= peak * (1 + tol)
        terminate = within_range or terminate
        if not terminate:
            newh, newq, newm = self.move_and_measure(range=90, steps=5, dry_run=False, a=maxX,
                                                b=maxY)
            h, q, m = data
            h = np.append(h, newh)
            q = np.append(q, newq)
            m = np.append(m, newm)
            data = h, q, m

            checkMaxH, checkMaxQ, checkMaxM = self.max_with_tolerance(*data, peak, tol)
            if checkMaxH is not None and checkMaxQ is not None:
                r0.move_to(maxX)
                r1.move_to(maxY)
                r0.wait_for_stop()
                r1.wait_for_stop()
                test_measurement = self.measure(measurements=1000)
                sleep(1)
                print((abs(test_measurement - peak) / peak))
                within_range = test_measurement >= peak * (1 - tol) and test_measurement <= peak * (1 + tol)
                if (within_range):
                    return X, Y, Z1, checkMaxH, checkMaxQ, checkMaxM
            return self.optimize(x0=None, data=data, terminate=True, tol=tol)
        else:
            return X, Y, Z1, maxX, maxY, maxZ
    def close(self):
        self.stage[0].close()
        self.stage[1].close()
        self.daq_task.close()
        print("Closed successfully")



class RotorExperiment(EnvExperiment):

    def run(self):
        rotor1 = RotatorFeedbackChannel(ch_name="Dev1/ai0", rotator_sn=["55000741", "55105674"], dry_run=False, plate_config = (qwp, hwp, arb_retarder))
        arb_angs = (30,45,23)
        q_angs = np.random.rand(5)*180
        h_angs = np.random.rand(5)*180
        output = rotor1.measure_sim(angles=(q_angs, h_angs, arb_angs), E = 2, background = 0)
        print(output[1])
        init_pos_1 = rotor1._pos(rotor_num=0)
        init_pos_2 = rotor1._pos(rotor_num=1)
        range_val = 180
        steps = 5
        phi_0 = None
        theta_0 = None
        eta_0 = None
        E_0 = None
        data = rotor1.move_and_measure(phi=phi_0, theta=theta_0, eta=eta_0, range=range_val,
                                steps=steps, dry_run=False, E=E_0, a=init_pos_1, b=init_pos_2)

        X, Y, Z1, maxX, maxY, maxZ = rotor1.optimize(data=data, alpha=init_pos_1, beta=init_pos_2,
                                              tol=0.05)
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






