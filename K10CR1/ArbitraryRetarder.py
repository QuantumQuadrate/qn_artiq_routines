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
    ar00,  ar01, ar10, ar11 = arb_retarder(args, piecewise=True)

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

todo: resturcture to allow scipy.minimize to use
"""
def plate_config_measure(configs):
    def generated_func(angles, input=np.array([1, 0]), E = 1, background = 0):
        output0, output1 = input.astype(complex)
        temp_output0 = output0
        temp_output1 = output1
        for c, a in zip(configs, angles):
            config00, config01, config10, config11 = c(a, piecewise = True)
            temp_output0 = config00 * output0 + config01 * output1
            temp_output1 = config10 * output0 + config11 * output1
            output0 = temp_output0
            output1 = temp_output1
        return np.sqrt(abs(output0**2))*E + background

    return generated_func


#Constrains that the max found by the minimize function must be greater or equal to the max measured value

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


