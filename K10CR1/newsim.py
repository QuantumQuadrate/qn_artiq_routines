import numpy as np
from scipy.optimize import minimize, curve_fit
import matplotlib.pyplot as plt
from time import sleep

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

def gen_phi(default = False):
    if not default:
        phi_x = np.random.rand(1) * np.pi/2
        phi_y = np.random.rand(1) * np.pi/2
    return phi_x, phi_y
def gen_state(default = False, phi_x = None, phi_y = None):
    REAL_STATE = np.array([1,0])

    if phi_x is None or phi_y is None:
        phi_x = np.random.rand(1) * np.pi/2
        phi_y = np.random.rand(1) * np.pi/2

    cx = np.cos(phi_x)
    sx = np.sin(phi_x)

    cy = np.cos(phi_y)
    sy = np.sin(phi_y)


    input_state = [cx - 1j*sx, cy - 1j*sy]
    if default:
        return REAL_STATE
    else:
        input_state /= np.linalg.norm(input_state)
        return input_state

def gen_secrets():
    phi_x, phi_y = gen_phi()
    return phi_x, phi_y, np.random.rand(1)*90
def move_and_measure(phi_x, phi_y, rand_axis, range, steps, a = 90, b = 0):
    q_ang = np.concatenate((( np.linspace(a, a+range, steps)), (np.linspace(a+range, a, steps ))))
    h_ang = np.concatenate((((np.linspace(b, b+range, steps)), (np.linspace(b+range, b+2*range, steps)))))
    measurements = []
    for q,h in zip(q_ang, h_ang):
            measurements.append(np.sum(measure(q_ang=q, h_ang=h, ax_trans= rand_axis, phi_x=phi_x, phi_y = phi_y)
                            ))

    return q_ang, h_ang, measurements

def measure(q_ang = 45, h_ang = 100, ax_trans = 75, phi_x = None, phi_y = None):
    qwp00, qwp01, qwp10, qwp11 = qwp(fast_axis_angle=q_ang, piecewise=True)
    hwp00, hwp01, hwp10, hwp11 = hwp(fast_axis_angle=h_ang, piecewise=True)
    lp00,  lp01, lp10, lp11 = lp(ax_trans=ax_trans,piecewise=True)
    input_state = gen_state(phi_x=phi_x, phi_y = phi_y)

    A_p = qwp00*input_state[0] + qwp01*input_state[1]
    B_p = qwp10*input_state[0] + qwp11*input_state[1]
    A_pp = hwp00*A_p + hwp01*B_p
    B_pp = hwp10*A_p + hwp11*B_p
    state_lp1 = lp00 * A_pp + lp01 * B_pp
    state_lp2 = lp10 * A_pp + lp11 * B_pp
    value = np.sqrt(abs(state_lp1**2 + state_lp2)**2)

    return value

def objective_func(x, args):
    print(x)
    ax_trans = x[0]
    phi_x = x[1]
    phi_y = x[2]
    q_ang = args[0]
    h_ang = args[1]
    measurements = args[2]
    error = 0
    #rint(args)
    for q, h, measurement in zip(q_ang, h_ang, measurements):
        error += np.sum(measure(q_ang=q, h_ang =h, ax_trans=ax_trans,phi_x = phi_x, phi_y = phi_y)-measurement)
    return error*1000

phi_x, phi_y, rand_axis = gen_secrets()
print(phi_x, phi_y, rand_axis)
q_ang, h_ang, measurements = move_and_measure(phi_x=phi_x, phi_y = phi_y, rand_axis=rand_axis, range = 45, steps = 1000)

initial_guess = np.random.rand(3)*np.pi -np.pi/2
#print(initial_guess)
args = [[],[],[]]
args[0].append(q_ang)
args[1].append(h_ang)
args[2].append(measurements)
#print(args)
bounds = [(0,np.pi), (0,np.pi), (0,180)]
result = minimize(objective_func, x0 = initial_guess, args=args, bounds=bounds, method='L-BFGS-B', tol = 1e6,
                  options={"maxit":1e20})
print(result)
x = result.x

phi_x1, phi_y1,ax_trans,  = x

predicted_state = gen_state(default=False, phi_x=phi_x, phi_y= phi_y)
print(predicted_state)
q = np.linspace(0,180, 100)
h = np.linspace(0,180, 100)
X, Y = np.meshgrid(q,h)

Z = measure(q_ang = X, h_ang = Y, phi_x=phi_x, phi_y=phi_y, ax_trans=rand_axis)
Z1 = measure(q_ang = X, h_ang = Y, phi_x=phi_x1, phi_y=phi_y1, ax_trans=ax_trans)
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.contour3D(X, Y, Z1, 50, cmap='binary')
ax.contour3D(X, Y, Z, 50,)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_title('3D contour')
plt.show()