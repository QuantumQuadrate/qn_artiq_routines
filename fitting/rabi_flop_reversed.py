import numpy as np
from scipy.optimize import minimize_scalar
from scipy.signal import lombscargle
from . import FitBase
"""Fit a typical Rabi flop time scan with a decaying cosine curve, including
initial dead time to account for AOM/... switching effects.

rabi_flop.py only supports positive-sign scans (starting at y=1, not y=0).
rabi_flop_reversed.py supports negative-sign scans (starting at y=0, not y=1)

For guessing the initial parameters, it is assumed that the scan range contains
the first minimum (i.e., the pi time), and is about 1/10 to 10 times in length.
"""


def parameter_initialiser(x, y, p):
    t_min = np.amin(x)
    t_range = np.amax(x) - t_min
    if t_range == 0.0:
        t_range = 1.0

    # Estimate frequency. Starting with a Lomb-Scargle periodogram (which
    # supports irregularly-spaced samples), we pick the strongest frequency
    # component which leads to a pi time larger than t_min.
    #
    # TODO: Could use better heuristics for frequency range based on minimum
    # distance between points -> aliasing.
    freq = np.pi / t_range
    freqs = np.linspace(0.1 * freq, 10 * freq, 2 * len(x))
    pgram = lombscargle(x, y, freqs, precenter=True)
    freq_order = np.argsort(-pgram)
    for f in freqs[freq_order]:
        t = 2 * np.pi / f
        if t / 2 > t_min:
            p["t_period"] = t
            break

    p["t_dead"] = 0.0

    p["y_lower"] = np.clip(2 * np.mean(y) - 1, 0, 1)

    # TODO: Estimate decay time constant using RMS amplitude from global mean
    # in first and last chunk.
    p["tau_decay"] = 1


def fitting_function(x, p):
    y_upper = 1.0
    shifted_t = (x - p["t_dead"])
    y = p["y_lower"] + (y_upper - p["y_lower"]) / 2 * (
        np.exp(-shifted_t / p["tau_decay"]) *
        np.cos(-np.pi + 2 * np.pi / p["t_period"] * shifted_t) + 1)
    return np.where(x < p["t_dead"], y_upper, y)


def derived_parameter_function(p, p_err):
    non_decaying_pi_time = p["t_dead"] + p["t_period"] / 2

    # Compute the point of maximum population transfer (minimum in y) which
    # will be slightly shifted towards zero in the face of non-zero tau_decay.
    fit = minimize_scalar(lambda t: fitting_function(t, p),
                          method="brent",
                          bracket=[0.9 * non_decaying_pi_time, non_decaying_pi_time])
    # if fit.success:
    #     p["t_pi"] = fit.x
    # else:
    #     p["t_pi"] = non_decaying_pi_time

    p["t_pi"] = non_decaying_pi_time

    # This is just a Gaussian error propagation guess.
    p_err["t_pi"] = np.sqrt(p_err["t_dead"]**2 + (p_err["t_period"] / 2)**2)
    return p, p_err


rabi_flop_reversed = FitBase.FitBase(
    ["t_period", "t_dead", "y_lower", "tau_decay"],
    fitting_function,
    parameter_initialiser=parameter_initialiser,
    derived_parameter_function=derived_parameter_function,
    parameter_bounds={
        "t_period": (0, np.inf),
        "t_dead": (0, np.inf),
        "y_lower": (0, 1),
        "tau_decay": (0, np.inf)
    })



# import numpy as np
# from . import FitBase
#
#
# def parameter_initialiser(x, y, p):
#     p["r"] = 0.9
#     p["C"] = 0.9
#     p["f"] = 10e3
#     p["x0"] = 0.0
#     p["tau"] = 1e-5
#
# def fitting_function(x, p):
#     y = (p["C"]/2)*np.cos(-np.pi + 2*np.pi*p["f"]*(x-p["x0"]))*np.exp(-(x-p["x0"])/abs(p["tau"])) + (p["r"] - p["C"]/2)
#
#     # y = np.where(x <= p['x0'], p['k1'] * (x - p['x0']) + p['d'],
#     #              p['k2'] * (x - p['x0']) + p['d'])
#     return y
#
#
# rabi_flop_reversed = FitBase.FitBase(['r', 'C', 'f', 'x0', 'tau'],
#                              fitting_function,
#                              parameter_initialiser=parameter_initialiser)
