import numpy as np
from scipy.optimize import minimize
from scipy.signal import argrelmax, argrelmin, argrelextrema
import matplotlib.pyplot as plt
from time import sleep
#from rotator_feedback import RotatorFeedbackChannel

rotor1 = 1
rotor2 = 2

def rotated(zero, one, two, three, theta):
    r00, r01, r10, r11 = np.cos(theta), np.sin(theta), np.cos(theta), np.cos(theta)
    p_00 = r00*zero + r10*one
    p_01 = r01*zero + r11*one
    p_10 = r00*two + r10*three
    p_11 = r01*two + r11*three

def qwp(fast_axis_angle = 90, piecewise = False):

    rad_angle = np.radians(fast_axis_angle)
    complex_factor = np.exp(complex(-1j * np.pi / 4))
    c = np.cos(rad_angle)
    s = np.sin(rad_angle)
    mat = np.array([[c**2 + 1j*s**2, (1-1j)*s*c],
                        [(1-1j)*s*c, s**2+1j*c**2]])
    matrix = mat.astype(complex)
    matrix *= complex_factor
    if not piecewise:
        return matrix
    else:
        return matrix[0][0], matrix[0][1], matrix[1][0], matrix[1][1]

def hwp(fast_axis_angle = 90, piecewise = False):
    rad_angle = np.radians(fast_axis_angle)
    complex_factor = np.exp(complex(-1j * np.pi / 2))
    c = np.cos(rad_angle)
    s = np.sin(rad_angle)
    mat = np.array([[c**2 - s**2, 2*s*c],
                        [2*s*c, s**2 - c**2]])
    matrix = mat.astype(complex)
    matrix *= complex_factor
    if not piecewise:
        return matrix
    else:
        return matrix[0][0], matrix[0][1], matrix[1][0], matrix[1][1]

def lp(ax_trans = 180, piecewise = False):
    rad_angle = np.radians(ax_trans)
    c = np.cos(rad_angle)
    s = np.sin(rad_angle)
    matrix = np.array([[c ** 2, c * s],
                       [c*s, s ** 2]])
    if not piecewise:
        return matrix
    else:
        return matrix[0][0], matrix[0][1], matrix[1][0], matrix[1][1]

def gen_phi(default = False):
    if not default:
        phi_x = np.random.rand(1) * 180
        phi_y = np.random.rand(1) * 180
    else:
        phi_x, phi_y = 0, 90
    return phi_x, phi_y
def gen_state(default = True, phi_x = None, phi_y = None, E = 1):
    REAL_STATE = np.array([[1,0]])

    phi_x = np.radians(phi_x)
    phi_y = np.radians(phi_y)

    if phi_x is None or phi_y is None:
        phi_x = float(np.random.rand(1)) * np.pi/2
        phi_y = float(np.random.rand(1)) * np.pi/2

    cx = np.cos(phi_x)
    sx = np.sin(phi_x)

    cy = np.cos(phi_y)
    sy = np.sin(phi_y)


    input_state = [cx - 1j*sx, cy - 1j*sy]
    if default:
        input_state /= np.linalg.norm(input_state)
        return [input_state[0]*E, 0]
    else:
        input_state /= np.linalg.norm(input_state)

        return [input_state[0]*E, 0]

def gen_secrets(default = True, E = None):
    phi_x, phi_y = gen_phi(default=default)
    axis = -1
    if E is None:
        E = np.random.rand()*2
    if default:
        axis = 45
        theta_h = 45
        theta_q = 45
    else:
        axis = np.random.rand(1)*180 - 90
        theta_q = np.random.rand(1) * 180 - 90
        theta_h = np.random.rand(1) * 180 - 90
    return phi_x, phi_y, axis, theta_q,theta_h, E

def move_and_measure(phi_x, phi_y, rand_axis, range, steps, a = 0, b = 0, theta_h = 0, theta_q = 0, dry_run = True,
                     r_feedback = None, E = None):
    q_ang = np.random.uniform(low = -range/2,high = range/2,size =steps) + a - theta_q
    h_ang = np.random.uniform(low = -range/2, high = range/2, size = steps) + b - theta_h
    q_ang.sort()
    h_ang.sort()
    measurements = []
    if dry_run or r_feedback is None:
        for q,h in zip(q_ang, h_ang):
            measurements.append(np.sum(measure(q_ang=q, h_ang=h, ax_trans= rand_axis, phi_x=phi_x, phi_y = phi_y, E = E)
                            ))
    else:
        for q, h in zip(q_ang, h_ang):
            r_feedback.stage[1].move_to(q)
            r_feedback.stage[0].move_to(h)
            r_feedback.stage[0].wait_for_stop()
            r_feedback.stage[1].wait_for_stop()
            measurement = r_feedback.measure()
            measurements.append(measurement)
            sleep(1)
            print(measurement)

    return q_ang, h_ang, measurements

def measure(q_ang = 45, h_ang = 100, ax_trans = 75, phi_x = None, phi_y = None, E = 1,  theta_h = 0, theta_q = 0):
    qwp00, qwp01, qwp10, qwp11 = qwp(fast_axis_angle=q_ang-theta_q, piecewise=True)
    hwp00, hwp01, hwp10, hwp11 = hwp(fast_axis_angle=h_ang-theta_h, piecewise=True)

    lp00,  lp01, lp10, lp11 = lp(ax_trans=ax_trans,piecewise=True)
    input_state = gen_state(phi_x=phi_x, phi_y = phi_y, E = E)

    A_p = qwp00*input_state[0] + qwp01*input_state[1]
    B_p = qwp10*input_state[0] + qwp11*input_state[1]
    A_pp = hwp00*A_p + hwp01*B_p
    B_pp = hwp10*A_p + hwp11*B_p
    state_lp1 = lp00 * A_pp + lp01 * B_pp
    state_lp2 = lp10 * A_pp + lp11 * B_pp
    value = np.sqrt(abs(state_lp1)**2 + abs(state_lp2)**2)

    return value

def objective_func(x, args):
    print(x)
    ax_trans = x[0]
    phi_x = x[1]
    phi_y = x[2]
    E = x[3]
    phase_q = x[4]
    phase_h = x[5]
    q_ang = args[0] - phase_q
    h_ang = args[1] - phase_h
    measurements = args[2]
    error = 0
    #print(args)
    for q, h, measurement in zip(q_ang, h_ang, measurements):

        error += (np.sum(measure(q_ang=q, h_ang =h, ax_trans=ax_trans,phi_x = phi_x, phi_y = phi_y, theta_q=phase_q,
                                 theta_h=phase_h, E=E)
                         -measurement))**2
    return error*1000
range_val = 90
steps = 20
phi_x, phi_y, rand_axis, theta_q, theta_h, E = gen_secrets(default=False, E = 2)
print(phi_x, phi_y, rand_axis)


#rotor_channel = RotatorFeedbackChannel(ch_name="Dev1/ai0", rotator_sn=["55105674", "55000741"], dry_run=False)
#rotor_channel.stage[0].stop()
#rotor_channel.stage[1].stop()
#rotor_channel.stage[0].move_to(0)
#rotor_channel.stage[1].move_to(0)
#rotor_channel.stage[0].wait_for_stop()
#rotor_channel.stage[1].wait_for_stop()

q_ang, h_ang, measurements = move_and_measure(phi_x=phi_x, phi_y = phi_y, rand_axis=rand_axis, range = range_val,
                                              steps = steps, r_feedback = None, dry_run=True, E = E,
                                              a=0,
                                              b= 0)

initial_guess = [np.random.rand()*180-90,np.random.rand()*180-90,np.random.rand()*180-90, max(measurements), 0, 0]


args = [[],[],[]]

args[0].append(q_ang)
args[1].append(h_ang)
args[2].append(measurements)
#print(args)

bounds = [(-90, 90), (-90, 90), (-90, 90),(max(measurements),np.inf), (0,180), (0,180)]

result = minimize(objective_func, x0 = initial_guess, args=args, bounds=bounds, method='trust-constr', tol = 1e-6,)

x = result.x
phi_x2, phi_y2, ax_trans, E_predict, theta_q2, theta_h2  = x

#print(result)
result2 = minimize(objective_func, x0 = x, args=args, bounds=bounds, method='trust-constr', tol = 1e-6,)
#print(result2)

actual_state = gen_state(default=False, phi_x=phi_x, phi_y=phi_y, E = E)
predicted_state = gen_state(phi_x = phi_x2, phi_y=phi_y2, E=E_predict)
print(predicted_state)
print(actual_state)
q = np.linspace(-range_val,range_val, 30)
h = np.linspace(-range_val,range_val, 30)
X, Y = np.meshgrid(q,h)

Z = measure(q_ang=X, h_ang=Y, phi_x=phi_x, phi_y=phi_y, ax_trans=rand_axis, theta_q = theta_q, theta_h = theta_h, E = E)
Z1 = measure(q_ang=X, h_ang=Y, phi_x=phi_x2, phi_y=phi_y2, ax_trans=ax_trans, theta_q = theta_q2, theta_h = theta_h2, E = E_predict)


#Z_scatter = measure(q_ang=q, h_ang=h, phi_x=phi_x, phi_y=phi_y, ax_trans=rand_axis)
#Z_1scatter = measure(q_ang=q, h_ang=h, phi_x=phi_x, phi_y=phi_y, ax_trans=ax_trans)

fig = plt.figure()
#ax = plt.axes(projection='3d')
ax = fig.add_subplot(projection ="3d")
ax.scatter(q_ang, h_ang, measurements, s = 30, marker="*")
ax.scatter(X, Y, Z, marker=".", s=50)
ax.scatter(X, Y, Z1, marker = ".", s=5)
c = argrelmax(Z1, order = 100)
maxX, maxY, maxZ = X[c], Y[c], Z1[c]
measurementMaxIndex = np.where(Z1 == (max(Z1[c])))
print(maxX[measurementMaxIndex[0]])
print((maxY[measurementMaxIndex[0]]))

#q_ang, h_ang, measurements = move_and_measure(phi_x=phi_x, phi_y = phi_y, rand_axis=rand_axis, range = range_val, steps = 10)

ax.scatter(maxX,maxY,maxZ, marker = "*", s = 20)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_title(' ')
plt.show()
