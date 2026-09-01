"""
example usage for the temp fitting code

This file is a reading example, not a runnable script: it references
placeholder data (t_steps_us, retention, errs, ...) that you supply in your
own analysis. The body is wrapped in a function so the ARTIQ repository
scan can import the file without executing it.
"""


def _example_usage(t_steps_us, retention, errs, temp_guess_uK,
                   baseline_retention_guess):
    # import the fitting functions
    import sys
    import numpy as np
    from matplotlib import pyplot as plt
    sys.path.append("C:\\..\\qn_artiq_routines")
    from fitting.run_modeling import get_release_recap_fit_result, release_recap_retention_at_t

    # define trap params
    trap_params = {'Tdepth': 2e-3, 'wx': 0.8e-6, 'lmda': 8.52e-7}

    # get the fit result for temp (K) and retention, and the retention points for the fit at each timestep.
    # pass in the time steps (us), retention, and a guess for the temp (uK) and baseline_retention
    popt, fit_retention = get_release_recap_fit_result(t_steps_us, retention, p0=[temp_guess_uK,baseline_retention_guess],
                                                       retention_at_t_kwargs=trap_params)

    # generate a temp curve with more timesteps since we typically don't have more than a few timesteps for the data
    hi_res_t_steps_us = np.linspace(0, 100, 100)
    fit_retention_hi_res = release_recap_retention_at_t(hi_res_t_steps_us, T=popt[0], base_retention=popt[1], **trap_params)

    fig, ax = plt.subplots()
    ax.scatter(t_steps_us, retention)
    ax.errorbar(t_steps_us, retention, errs, ls='none')
    ax.plot(hi_res_t_steps_us, fit_retention_hi_res, label=f'fit: T={popt[0]*1e6:.2f} uK, r={popt[1]:.2f}',
            linestyle='dashdot')
