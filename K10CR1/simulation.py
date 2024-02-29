import numpy as np
from scipy.optimize import minimize, curve_fit
from scipy.signal import argrelmax
import matplotlib.pyplot as plt
from time import sleep

sin = np.sin
cos = np.cos
exp = np.exp
pi = np.pi
rand = np.random.rand


def exi(exponentiation: float):
    return exp(1j * exponentiation)
def rotation(theta):
    R00 = cos(theta)
    R01 = sin(theta)
    R10 = -1 * sin(theta)
    R11 = cos(theta)
    return R00, R01, R10, R11

def arb_retarder(phi: float, theta: float, eta: float):

    c = cos(theta)
    s = sin(theta)
    gp = exi(-eta/2)

    off_diag_factor = (1 - exi(eta)) * c * s
    J00 = c**2 + exi(eta)*(s**2)
    J01 = off_diag_factor*exi(-1*phi)
    J10 = off_diag_factor*exi(phi)
    J11 = exi(eta) *c ** 2 +  (s ** 2)

    J00 *= gp
    J01 *= gp
    J10 *= gp
    J11 *= gp

    return J00, J01, J10, J11


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

def measure(q_ang = 45, h_ang = 100, theta = 3*pi/4, phi = pi/3, eta = pi/4, E = 1,  theta_h = 0, theta_q = 0):
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

    return value


def move_and_measure(theta, eta, phi, range, steps, a = 0, b = 0, theta_h = 0, theta_q = 0, dry_run = True,
                     r_feedback = None, E = None):
    q_ang = np.random.uniform(low = -range/2,high = range/2,size =steps) + a - theta_q
    h_ang = np.random.uniform(low = -range/2, high = range/2, size = steps) + b - theta_h
    q_ang.sort()
    h_ang.sort()
    measurements = []
    if dry_run or r_feedback is None:
        for q,h in zip(q_ang, h_ang):
            measurements.append(np.sum(measure(q_ang=q, h_ang=h, theta = theta, phi=phi,
                                               eta = eta, E = E)))
    """else:
        for q, h in zip(q_ang, h_ang):
            r_feedback.stage[1].move_to(q)
            r_feedback.stage[0].move_to(h)
            r_feedback.stage[0].wait_for_stop()
            r_feedback.stage[1].wait_for_stop()
            measurement = r_feedback.measure()
            measurements.append(measurement)
            sleep(1)
            print(measurement)"""

    return q_ang, h_ang, measurements


def objective_func(x, args):
    theta = x[0]
    eta = x[1]
    phi = x[2]
    E = x[3]
    phase_q = 0
    phase_h = 0
    q_ang = args[0]
    h_ang = args[1]
    measurements = args[2]
    error = 0

    for q, h, measurement in zip(q_ang, h_ang, measurements):

        error += (np.sum(measure(q_ang=q, h_ang =h, theta=theta, phi = phi, eta = eta, theta_q=phase_q,
                                 theta_h=phase_h, E=E)
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

range_val = 90
steps = 30
theta_0, eta_0, phi_0, E_0 = gen_secrets(default=False)

q_ang, h_ang, measurements = move_and_measure(phi=phi_0, theta = theta_0, eta = eta_0, range = range_val,
                                              steps = steps, r_feedback = None, dry_run=True, E = E_0, a=0, b= 0)

initial_guess = [rand()*pi, rand()*pi, rand()*pi, rand()*2]

args = [[],[],[]]

args[0].append(q_ang)
args[1].append(h_ang)
args[2].append(measurements)
#print(args)

bounds = [(0, pi), (0, pi), (0, pi), (max(measurements), np.inf)]
result = minimize(objective_func, x0 = initial_guess, args=args, bounds=bounds, method='trust-constr', tol = 1e-15,)
x = result.x
theta, eta, phi, E = x

q = np.linspace(-range_val,range_val, 30)
h = np.linspace(-range_val,range_val, 30)
X, Y = np.meshgrid(q,h)
print(E-E_0)
print(theta-theta_0)
print(eta-eta_0)
print(phi-phi_0)


Z = measure(q_ang=X, h_ang=Y, theta = theta_0, eta=eta_0, phi=phi_0, theta_q = 0, theta_h = 0, E = E_0)
Z1 = measure(q_ang=X, h_ang=Y, theta = theta, eta=eta, phi=phi, theta_q = 0, theta_h = 0, E = E)

fig = plt.figure()
#ax = plt.axes(projection='3d')
ax = fig.add_subplot(projection ="3d")
#ax.scatter(q_ang, h_ang, measurements, s = 30, marker="*")
ax.scatter(X, Y, Z, marker=".", s=3)
ax.scatter(X, Y, Z1, marker = ".", s=5)
c = argrelmax(Z1, order = 100)
#maxX, maxY, maxZ = X[c], Y[c], Z1[c]
#measurementMaxIndex = np.where(Z1 == (max(Z1[c])))
#print(maxX[measurementMaxIndex[0]])
#print((maxY[measurementMaxIndex[0]]))

#q_ang, h_ang, measurements = move_and_measure(phi_x=phi_x, phi_y = phi_y, rand_axis=rand_axis, range = range_val, steps = 10)

#.scatter(maxX,maxY,maxZ, marker = "*", s = 20)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_title(' ')
plt.show()