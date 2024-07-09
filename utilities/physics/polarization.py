"""
vectors, matrices, and functions for use in polarization calculations and visualization

Adapted from https://github.com/prhuft/rubidium/blob/master/optics/polarization_rotation_and_correction.ipynb,
which contains several usage examples.
"""

import numpy as np
from numpy import sin, cos, exp, pi, sqrt, trace, dot, vdot
import matplotlib.pyplot as plt

# states defined in the basis {1,0} = {H,V}
V = np.array([0, 1])
H = np.array([1, 0])
D = np.array([1, 1]) / sqrt(2)
A = np.array([1, -1]) / sqrt(2)
R = np.array([1, -1j]) / sqrt(2)
L = np.array([1, 1j]) / sqrt(2)

# in the H,V basis
s0 = np.array([[1, 0], [0, 1]])
s1 = np.array([[0, 1], [1, 0]])
s2 = np.array([[0, -1j], [1j, 0]])
s3 = np.array([[1, 0], [0, -1]])

# we want to switch to the R,L basis, which is more common for drawing on the Poincaré sphere
projR = np.outer(R, np.conj(R))
evals, evecs = np.linalg.eig(projR)
U = evecs.transpose()  # doesn't do anything in this case b.c. the matrix is symmetric about the diagonal
Uinv = np.linalg.inv(U)
assert np.isclose(trace(Uinv.dot(projR).dot(U)), np.sum(Uinv.dot(projR).dot(U))), "oops, matrix isn't diagonal"

# in the R,L basis
s1 = Uinv.dot(s1).dot(U)
s2 = Uinv.dot(s2).dot(U)
s3 = Uinv.dot(s3).dot(U)


def poincare_sphere(ax=None):
    """
    plots the Poincare sphere and returns the associated Axes object ax.

    The Z axis is vertical, as is standard, with circular polarization at the Z poles, H/V along X, and D/A along Y.

    ax: an Axes object. if not provided, one will be created
    """

    # Create a sphere
    phi, theta = np.mgrid[0:2 * np.pi:100j, 0:np.pi:50j]
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)

    if ax is None:
        fig = plt.figure(figsize=(5, 5), dpi=100)
        ax = fig.add_subplot(111, projection='3d')

    # Plotting the Poincaré sphere
    ax.plot_surface(x, y, z, color='pink', alpha=0.3, rstride=1, cstride=1, linewidth=0)

    # Labeling the poles
    ax.text(0, 0, 1.1, 'R', color='black', fontsize=14, ha='center')
    ax.text(0, 0, -1.1, 'L', color='black', fontsize=14, ha='center')

    ax.text(1.1, 0, 0, 'V', color='black', fontsize=14, ha='center')
    ax.text(-1.1, 0, 0, 'H', color='black', fontsize=14, ha='center')

    ax.text(0, 1.1, 0, 'D', color='black', fontsize=14, ha='center')
    ax.text(0, -1.1, 0, 'A', color='black', fontsize=14, ha='center')

    # Generate points for the equatorial circle
    theta = np.linspace(0, 2 * np.pi, 100)
    x_eq = np.cos(theta)
    y_eq = np.sin(theta)
    z_eq = np.zeros_like(theta)

    # Plot the equatorial circle
    ax.plot(x_eq, y_eq, z_eq, 'k--')  # 'k--' specifies a black dashed line

    axis_pts = np.linspace(-1, 1, 100)
    ax.plot(axis_pts, np.zeros_like(axis_pts), np.zeros_like(axis_pts), 'k--')
    ax.plot(np.zeros_like(axis_pts), axis_pts, np.zeros_like(axis_pts), 'k--')
    ax.plot(np.zeros_like(axis_pts), np.zeros_like(axis_pts), axis_pts, 'k--')

    # Setting the aspect ratio
    ax.set_box_aspect([1, 1, 1])

    # Hide the axes
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])

    # Transparent spines
    ax.xaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))

    # Transparent panes
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

    ax.grid(False)
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.set_aspect('equal')

    return ax


def get_stokes_params(jones_vec):
    """
    compute the four Stokes parameters from a 2-element Jones vector
    """
    rho = np.outer(np.conj(jones_vec), jones_vec)
    S0 = trace(dot(s0, rho))
    S1 = trace(dot(s1, rho))
    S2 = trace(dot(s2, rho))
    S3 = trace(dot(s3, rho))

    return np.real(S0), np.real(S1), np.real(S2), np.real(S3)


def draw_stokes_vector(ax, jones_vec, color='k', label=''):
    """
    plot the 2-element complex jones vector on the supplied ax

    ax is presumed to be a 3D projection.

    Note that the Stokes parameters S3, S2, and S1 used here are with respect to a basis {1,0} = {R,L}, which are
    are typically taken to be at the poles of, e.g. the Bloch sphere, or equivalently, along the S3 axis of the Poincare sphere.
    However, unlike the typical convention, {H,V} is taken to be along S2 (Y) and {D,A} along S1 (X), which is unlike
    the typical convention. I did this because I found that those were the components of the respective states after transforming
    the Pauli matrices to the R,L basis.

    return S0, S1, S2, S3
    """

    S0, S1, S2, S3 = get_stokes_params(jones_vec)

    ax.text(1.5 * S2, 1.5 * S1, 1.5 * S3, label, color=color, fontsize=14, ha='center')
    ax.quiver(0, 0, 0, S2, S1, S3, color=color, arrow_length_ratio=0.1)


def QWP(theta):
    """
    Quarter waveplate matrix.
    :param theta: angle of fast axis wrt vertical in radians
    :return: the 2x2 array qwp_at_theta
    """
    qwp_at_theta = np.array([[(0.5 - 0.5 * 1j) * (1j + cos(2 * theta)), (1j - 1) * cos(theta) * sin(theta)],
                             [(1j - 1) * cos(theta) * sin(theta), 1j * cos(theta) ** 2 + sin(theta) ** 2]])
    return qwp_at_theta


def HWP(theta):
    """
    Half waveplate matrix.
    :param theta: angle of fast axis wrt vertical in radians
    :return: the 2x2 array hwp_at_theta
    """
    hwp_at_theta = np.array([[cos(2 * theta), -2 * cos(theta) * sin(theta)],
                             [-2 * cos(theta) * sin(theta), -cos(2 * theta)]])
    return hwp_at_theta


def AWP(theta, phi, eta):
    """
    Arbitrary retarder unitary matrix parametrized by theta, phi, eta

    :param theta:
    :param phi:
    :param eta:
    :return: 2x2 array awp_at_theta_phi_eta
    """
    awp_at_theta_phi_eta = exp(-1j * eta / 2) * np.array([[cos(theta) ** 2 + exp(1j * eta) * sin(theta) ** 2,
                                                           (1 - exp(1j * eta)) * exp(-1j * phi) * cos(theta) * sin(
                                                               theta)],
                                                          [(1 - exp(1j * eta)) * exp(1j * phi) * cos(theta) * sin(
                                                              theta),
                                                           exp(1j * eta) * cos(theta) ** 2 + sin(theta) ** 2]])
    return awp_at_theta_phi_eta
