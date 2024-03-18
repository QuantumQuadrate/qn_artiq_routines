import numpy as np
from scipy.optimize import minimize, curve_fit
from scipy.signal import argrelmax
import matplotlib.pyplot as plt
from time import sleep
#from rotator_feedback import RotatorFeedbackChannel


sin = np.sin
cos = np.cos
exp = np.exp
pi = np.pi
rand = np.random.rand


def exi(exponentiation: complex):
    return exp(1j * exponentiation)
def rotation(theta):
    R00 = cos(theta)
    R01 = sin(theta)
    R10 = -1 * sin(theta)
    R11 = cos(theta)
    return R00, R01, R10, R11

def arb_retarder(phi: complex, theta: complex, eta: complex):
    c = cos(theta)
    s = sin(theta)
    gp = exi(-eta/2)
    off_diag_factor = (1 - exi(eta)) * c * s
    mat = np.array([[0,0],[0,0]]).astype(complex)
    mat[0][0] = c**2 + exi(eta)*(s**2)
    mat[0][1] = off_diag_factor*exi(-1*phi)
    mat[1][0] = off_diag_factor*exi(phi)
    mat[1][1] = exi(eta) *c ** 2 +  (s ** 2)

    mat[0][0] *= gp
    mat[0][1] *= gp
    mat[1][0] *= gp
    mat[1][1] *= gp

    return mat[0][0], mat[0][1], mat[1][0], mat[1][1]


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

def measure(q_ang = 45, h_ang = 100, theta = 3*pi/4, phi = pi/3, eta = pi/4, E = 1,  theta_h = 0, theta_q = 0, a = 0):
    qwp00, qwp01, qwp10, qwp11 = qwp(fast_axis_angle=q_ang-theta_q, piecewise=True)
    hwp00, hwp01, hwp10, hwp11 = hwp(fast_axis_angle=h_ang-theta_h, piecewise=True)

    ar00,  ar01, ar10, ar11 = arb_retarder(theta=theta, phi=phi, eta=eta)

    input_state = [E, 0]

    A_p = qwp00*input_state[0] + qwp01*input_state[1]
    B_p = qwp10*input_state[0] + qwp11*input_state[1]
    A_pp = hwp00*A_p + hwp01*B_p
    B_pp = hwp10*A_p + hwp11*B_p
    state_lp1 = ar00 * A_pp + ar01 * B_pp
    state_lp2 = ar10 * A_pp + ar11 * B_pp
    value = np.sqrt(abs(state_lp1)**2)

    return value + a


def move_and_measure(theta, eta, phi, range, steps, a = 0, b = 0, theta_h = 0, theta_q = 0, dry_run = True,
                     r_feedback = None, E = None):
    q_ang = np.random.uniform(low = -range/2,high = range/2,size =steps) + a - theta_q
    h_ang = np.random.uniform(low = -range/2, high = range/2, size = steps) + b - theta_h
    #q_ang.sort()
    #h_ang.sort()
    measurements = []
    if dry_run or r_feedback is None:
        for q,h in zip(q_ang, h_ang):
            measurements.append(np.sum(measure(q_ang=q, h_ang=h, theta = theta, phi=phi,
                                               eta = eta, E = E)))
    else:
        for q, h in zip(q_ang, h_ang):
            r_feedback.stage[1].move_to(q)
            r_feedback.stage[0].move_to(h)


            r_feedback.stage[1].wait_for_stop()
            r_feedback.stage[0].wait_for_stop()

            measurement = r_feedback.measure()
            measurements.append(measurement)
            sleep(1)
            print(measurement)

    return h_ang, q_ang, measurements

def cons1(x, args):
    return (x[3] - x[6] - args)
def cons2(x, h_arg,q_arg, m_arg):
    for h, q, m in zip(h_arg, q_arg, m_arg):
        error = 0
        return 0


def objective_func(x, hdata, qdata, mdata):
    error = 0

    for h, q, measurement in zip(hdata, qdata, mdata):
        print(q)
        error += (np.sum(measure(q_ang=q, h_ang = h, theta = x[0], eta = x[1], phi = [2], E = [3], theta_q=x[4], theta_h = x[5], a = x[6])
                         -measurement)*1000)**2
    return error*1000

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




rotor_channel = None #RotatorFeedbackChannel(ch_name="Dev1/ai0", rotator_sn=["55105674", "55000741"], dry_run=False)
#rotor_channel.stage[0].stop()
#rotor_channel.stage[1].stop()
#rotor_channel.stage[0].move_to(0)
#rotor_channel.stage[1].move_to(0)
#rotor_channel.stage[0].wait_for_stop()
#rotor_channel.stage[1].wait_for_stop()


range_val = 180
steps = 10
theta_0, eta_0, phi_0, E_0 = gen_secrets(default=False)

data = move_and_measure(phi=phi_0, theta = theta_0, eta = eta_0, range = range_val,
                                              steps = steps, r_feedback = rotor_channel, dry_run=True, E = E_0, a=0, b= 0)

pot_max = max(data[2])
pot_min = min(data[2])
x0 = [rand()*pi, rand()*pi, rand()*pi, pot_max*1.2, 0, 0, 0]


#print(args)
cons = ({'type': 'ineq',
       'fun': cons1,
       'args': (pot_max,)
       })

bounds = [(0, pi), (0, pi), (0, pi), (pot_max-pot_min, np.inf), (0, pi), (0, pi), (-0.5, pot_min)]
result = minimize(objective_func, x0 = x0, args=data, bounds=bounds, method='trust-constr', tol = 1e-15, constraints = cons)
x = result.x
theta, eta, phi, E, p_q, p_h, a = x

q = np.linspace(-range_val/2,range_val/2, 30)
h = np.linspace(-range_val/2,range_val/2, 30)
X, Y = np.meshgrid(q,h)
print(E-E_0)
print(theta-theta_0)
print(eta-eta_0)
print(phi-phi_0)


Z = measure(q_ang=X, h_ang=Y, theta = theta_0, eta=eta_0, phi=phi_0, theta_q = 0, theta_h = 0, E = E_0)
Z1 = measure(q_ang=X, h_ang=Y, theta = theta, eta=eta, phi=phi, theta_q = p_q, theta_h = p_h, E = E, a = a)

fig = plt.figure()
#ax = plt.axes(projection='3d')
ax = fig.add_subplot(projection ="3d")
#ax.scatter(q_ang, h_ang, measurements, s = 30, marker="*")
#ax.scatter(X, Y, Z, marker=".", s=3)
ax.scatter(X, Y, Z1, marker = ".", s=5)
ax.scatter(data[0], data[1], data[2])
c = argrelmax(Z1, order = 100)
maxX, maxY, maxZ = X[c], Y[c], Z1[c]
measurementMaxIndex = np.where(Z1 == (max(Z1[c])))
maxVals = np.sum(maxX[measurementMaxIndex[0]]), np.sum(maxY[measurementMaxIndex[0]])
#rotor_channel.move_to(degrees=maxVals[0], rotor_num=0)
#rotor_channel.move_to(degrees=maxVals[1], rotor_num=1)
#rotor_channel.wait_stop()
print(maxVals)

#q_ang, h_ang, measurements = move_and_measure(phi_x=phi_x, phi_y = phi_y, rand_axis=rand_axis, range = range_val, steps = 10)

#.scatter(maxX,maxY,maxZ, marker = "*", s = 20)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_title(' ')
plt.show()
rotor_channel.close()