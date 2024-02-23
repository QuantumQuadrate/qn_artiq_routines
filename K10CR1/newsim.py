import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
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
def gen_state(default = False):
    REAL_STATE = [[1/np.sqrt(2)],[1/np.sqrt(2)]]

    if default:
        return REAL_STATE
    else:
        input_state = np.random.rand(2) + 1j * np.random.rand(2)
        input_state /= np.linalg.norm(input_state)
        return input_state

def measure(ax_trans = 45, q_ang = 45, h_ang = 100, input_state = gen_state(default=True)):
    qwp00, qwp01, qwp10, qwp11 = qwp(fast_axis_angle=q_ang, piecewise=True)
    hwp00, hwp01, hwp10, hwp11 = hwp(fast_axis_angle=h_ang, piecewise=True)
    lp00,  lp01, lp10, lp11 = lp(ax_trans=ax_trans,piecewise=True)
    A_p = qwp00*input_state[0][0] + qwp01*input_state[1][0]
    B_p = qwp10*input_state[0][0] + qwp11*input_state[1][0]
    A_pp = hwp00*A_p + hwp01*B_p
    B_pp = hwp10*A_p + hwp11*B_p
    state_lp_measured = lp00 * A_pp + lp01 * B_pp

    value = abs(state_lp_measured**2)

    return value

q = np.linspace(0,360, 30)
h = np.linspace(0,360, 30)
X, Y = np.meshgrid(q,h)
Z = measure(q_ang = X, h_ang = Y, input_state=gen_state(default=False))
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.contour3D(X, Y, Z, 50, cmap='binary')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_title('3D contour')
plt.show()