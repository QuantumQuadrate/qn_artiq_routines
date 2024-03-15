import numpy as np
from scipy.optimize import minimize, curve_fit
from scipy.signal import argrelmax
import matplotlib.pyplot as plt
from time import sleep
from rotator_feedback import RotatorFeedbackChannel

sin = np.sin
cos = np.cos
exp = np.exp
pi = np.pi
rand = np.random.rand
###SOURCE: https://en.wikipedia.org/wiki/Jones_calculus


def exi(exponentiation: float):
    return exp(1j * exponentiation)


def rotation(theta):
    R00 = cos(theta)
    R01 = sin(theta)
    R10 = -1 * sin(theta)
    R11 = cos(theta)
    return R00, R01, R10, R11

"""
Function that returns the values for the Jones Matrix of an arbitrary retarding material at 
fast axis angle for x-axis being: theta
diff between fast and slow axis : eta
phi: the circularity
"""
def arb_retarder(args, piecewise = False):
    phi, theta, eta = args
    c = cos(theta)
    s = sin(theta)
    gp = exi(-eta/2)

    off_diag_factor = (1 - exi(eta)) * c * s
    J00 = c**2 + exi(eta)*(s**2)
    J01 = off_diag_factor*exi(-1*phi)
    J10 = off_diag_factor*exi(phi)
    J11 = exi(eta) *c ** 2 + (s ** 2)

    J00 *= gp
    J01 *= gp
    J10 *= gp
    J11 *= gp
    mat = [[J00, J01], [J10, J11]]
    if piecewise:
        return J00, J01, J10, J11
    else:
        return mat


"""
Function that returns the values for the Jones Matrix of the QWP at fast axis angle being fast_axis_angle
piecewise == True, returns the matrix as 4 seperate values,
which is used for computing alongside scipy.minimize, which
inputs multiple values at once, meaning you are required to do individual computation of 
matrix elements to preserve 
the dimensionality
"""
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


"""
Function that returns the values for the Jones Matrix of the HWP at fast axis angle being fast_axis_angle
piecewise == True, returns the matrix as 4 seperate values,
which is used for computing alongside scipy.minimize, which
inputs multiple values at once, meaning you are required to do individual computation of 
matrix elements to preserve 
the dimensionality 
"""
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

"""
Function simulating a measurement for light passing through a QWP at angle q_ang, then through HWP at angle h_ang, then 
through a arbitrarily retarding material, then measured along the x-axis 
"""
def measure(q_ang = 45, h_ang = 100, theta = 3*pi/4, phi = pi/3, eta = pi/4, E = 1,  theta_h = 0, theta_q = 0, a = 0):
    qwp00, qwp01, qwp10, qwp11 = qwp(fast_axis_angle=q_ang-theta_q, piecewise=True)
    hwp00, hwp01, hwp10, hwp11 = hwp(fast_axis_angle=h_ang-theta_h, piecewise=True)

    args  = theta, phi, eta
    ar00,  ar01, ar10, ar11 = arb_retarder(args)

    input_state = [1, 0]

    A_p = qwp00*input_state[0] + qwp01*input_state[1]
    B_p = qwp10*input_state[0] + qwp11*input_state[1]
    A_pp = hwp00*A_p + hwp01*B_p
    B_pp = hwp10*A_p + hwp11*B_p
    state_lp1 = ar00 * A_pp + ar01 * B_pp
    state_lp2 = ar10 * A_pp + ar11 * B_pp
    value = np.sqrt(abs(state_lp1)**2)

    return E*value + a

"""
#params:
#configs == 2x2 array of functions corresponding to a certain plate/retarder
#angles == arrays of inputs for the waveplates dictating their position
#input == 2x1 vector describing an EM wave in form [E0 exp(phi_x), E0 exp(phi_y)]
"""
def plate_config_measure(configs):
    def generated_func(angles, input=[1, 0]):
        output = input
        temp_output = output
        for c, a in zip(configs, angles):
            config = c(a)
            temp_output[0] = config[0][0] * output[0] + config[0][1] * output[1]
            temp_output[1] = config[1][0] * output[0] + config[1][1] * output[1]
            output = temp_output
        return output

    return generated_func
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

#Constrains that the max found by the minimize function must be greater or equal to the max measured value
def constraint(x, args):
    return (x[3] + x[6] - args)

configs = [hwp, qwp, arb_retarder]
arb  = [1.3,4.5, 3.2]
angles = [0, 45, arb]
input = [1,0]

gen_func = plate_config_measure(configs)

print(gen_func([0,45]))