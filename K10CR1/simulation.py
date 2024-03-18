import numpy as np
from scipy.optimize import minimize, curve_fit
from scipy.signal import argrelmax
import matplotlib.pyplot as plt
from time import sleep
from rotator_feedback import RotatorFeedbackChannel
import math
from ArbitraryRetarder import *


def q_h_gen(range, steps, center_x, center_y, random = False):

    if not random:
        range /= 2
        rads = np.arange(0, 2*np.pi, (2*np.pi)/steps)
        x = range * ( np.cos(rads) ) - center_x
        y = range * ( np.sin(rads) ) - center_y
        x1 = range/2 * ( np.cos(rads) ) - center_x
        y1 = range/2 * ( np.sin(rads) ) - center_y
        x2 = range / 4 * (np.cos(rads)) - center_x
        y2 = range / 4 * (np.sin(rads)) - center_y
        x0 = np.append(x, x1)
        x0 = np.append(x0, x2)
        y0 = np.append(y, y1)
        y0 = np.append(y0, y2)
    if random is True:
        x0, x1, x2 = rand(steps)*range, rand(steps)*range/2, rand(steps)*range/4
        x0 = np.append(x0,x1)
        x0 = np.append(x0, x2)
        y0, y1, y2 = rand(steps)*range, rand(steps)*range/2, rand(steps)*range/4
        y0 = np.append(y0, y1)
        y0 = np.append(y0, y2)

    return x0,y0
def move_and_measure(theta = None, eta = None, phi = None, range = None, steps = None, a = 0, b = 0, theta_h = 0, theta_q = 0, dry_run = True,
                     r_feedback = None, E = None):
    h_ang, q_ang = q_h_gen(steps=steps, range=range, center_x = -a, center_y= -b, random=True)
    q_ang -=  theta_q
    h_ang -=  theta_h
    ang1 = []
    ang2 = []

    measure = r_feedback.measure
    measurements = []
    if dry_run or r_feedback is None:
        for h,q in zip(q_ang, h_ang):
            measurements.append(np.sum(measure(q_ang=q, h_ang=h, theta = theta, phi=phi,
                                               eta = eta, E = E)))
    else:
        for a1, a2 in zip(q_ang, h_ang):
            r_feedback.stage[0].move_to(a1)
            r_feedback.stage[1].move_to(a2)
            i = 0
            sleep(0.3)
            while r_feedback.is_moving():
                mq_ang = np.append(mq_ang, r_feedback._pos(rotor_num=1))
                mh_ang = np.append(mh_ang, r_feedback._pos(rotor_num=0))
                measurement = measure(measurements=1)
                measurements = np.append(measurements, measurement)
                #print("intermediate measurement")
                #print(q_ang)
                i+=1
                if not i%10:
                    print("still reading...")
                sleep(0.3)
            r_feedback.stage[1].wait_for_stop()
            r_feedback.stage[0].wait_for_stop()


    return ang1, ang2, measurements


def measurement_constraint(x, h_args, q_args, m_args):
    error = 0
    for h, q, m in zip(h_args, q_args, m_args):
        error += m - measure(h, q ,*x)
    return error
def measurement_constraint2(x,h,q,m):
    return (measure(h,q,*x)-m)
#function to be minimized by scipy.minimize, emulates the behavior of some arbitrary light retarder
def max_constraint(x):
    return 0
def objective_func(x, hdata, qdata, mdata):
    theta = x[0]
    eta = x[1]
    phi = x[2]
    E = x[3]
    phase_q = x[4]
    phase_h = x[5]
    a = x[6]
    h_ang = hdata
    q_ang = qdata
    measurements = mdata
    error = 0

    for h, q, measurement in zip(h_ang, q_ang, measurements):

        error += (np.sum(measure(q_ang=q, h_ang =h, theta=theta, phi = phi, eta = eta, theta_q=phase_q,
                                 theta_h=phase_h, E=E, a = a)
                         -measurement))**2
    return error

def gen_secrets(default = False):
    if default:
        theta = pi/3
        eta = pi
        phi = 2*pi/2
        E = 1
        return theta, eta, phi, E
    else:
        theta = rand()*pi
        eta = rand()*pi
        phi = rand()*pi
        E = rand()*2

        return theta, eta, phi, E

def max_with_tolerance(h,q,m,peak,tol):
    m_indices = [i < (peak*(1+tol)) and i > (peak*(1-tol)) for i in m]
    m_indices = [j == max(m) and t for j, t in zip(m,m_indices)]
    for i, j, m in zip(h, q, m):
        if i != -361 and j != -361:
            return i, j, m
    return None, None, None
def optimize(m_func, r0, r1, data = None, x0 = None, bounds = None, rotor_channel = None, cons= None, range_val = 180, terminate = False, alpha = 0, beta = 0, tol = 0.2):
    cons = ({'type': 'ineq',
             'fun': constraint,
             'args': (max(data[2]),)
             },
            {'type': 'eq',
             'fun': max_constraint,
             },)


    method = 'trust-constr'
    trust_constr_opts = {'disp': True, 'barrier_tol':1e-8, 'xtol': 1e-12,
                         'initial_constr_penalty':1, 'maxiter':5e3,}
    if data is not None:
        pot_max = max(data[2])
        pot_min = min(data[2])
        bounds = [(0, pi), (0, pi), (0, pi), (pot_max - pot_min, np.inf), (0, pi), (0, pi),
                  (-0.5, pot_max)]
    else:
        pot_max = rand()*3
    if x0 is None:
        x0 = [rand() * pi, rand() * pi, rand() * pi, pot_max-pot_min, 0, 0, 0]
        sleep(1)
    result = minimize(objective_func, x0=x0, args = data, bounds=bounds, method=method,
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


    r0.move_to(maxX)
    r1.move_to(maxY)
    r0.wait_for_stop()
    r1.wait_for_stop()
    test_measurement = m_func(measurements=1000)
    sleep(1)
    print( (abs(test_measurement - peak) / peak ) )
    within_range = test_measurement >= peak*(1-tol) and test_measurement <= peak*(1+tol)
    terminate = within_range or terminate
    if not terminate:
        newh, newq, newm = move_and_measure(range=90, steps=5, r_feedback=rotor_channel, dry_run=False, a=maxX, b=maxY)
        h, q, m = data
        h = np.append(h, newh)
        q = np.append(q, newq)
        m = np.append(m, newm)
        data = h, q, m

        checkMaxH, checkMaxQ, checkMaxM = max_with_tolerance(*data, peak, tol)
        if checkMaxH is not None and checkMaxQ is not None:
            r0.move_to(maxX)
            r1.move_to(maxY)
            r0.wait_for_stop()
            r1.wait_for_stop()
            test_measurement = m_func(measurements=1000)
            sleep(1)
            print((abs(test_measurement - peak) / peak))
            within_range = test_measurement >= peak * (1 - tol) and test_measurement <= peak*(1 + tol)
            if(within_range):
                return X, Y, Z1, checkMaxH, checkMaxQ, checkMaxM
        return optimize(m_func=m_func, r0 = r0, r1=r1, x0 = None, data=data, rotor_channel=rotor_channel, terminate=True, tol = tol)
    else:
        return X, Y, Z1, maxX, maxY, maxZ


rotor_channel = RotatorFeedbackChannel(ch_name="Dev1/ai0", rotator_sn=["55000741", "55105674" ], dry_run=False, )

range_val = 180
steps = 4
theta_0, eta_0, phi_0, E_0 = gen_secrets(default=False)
