from artiq.experiment import *
from artiq.language import us, ns, MHz
import logging

import numpy as np
import unittest

import sys, os

cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")

from fitting_oxford.sinusoid import sinusoid

from utilities.BaseExperiment import BaseExperiment
from utilities.conversions import dB_to_V_kernel as dB_to_V
from subroutines.k10cr1_functions import *

logger = logging.getLogger(__name__)

class SinusoidTest(EnvExperiment):
    def build(self):
        self.base = BaseExperiment(experiment=self)
        self.base.build()
        self.base.set_datasets_from_gui_args()


    def prepare(self):
        self.base.prepare()


    def run(self):
        # p, p_err = self.test_random_data()
        # self.test_random_data_fixed_phi()
        self.test_pi_pulse()
        self.test_pi_pulse_with_t_dead()
        self.test_delayed_split_data_with_t_dead()


    def test_random_data(self):
        n_sample = 30
        t_max = 1e-1
        amp = 1.5
        rel_noise = 0.1
        offset = 40
        phi = np.pi / 2
        t_dead = 0.0
        omega = np.pi / t_max * 8.2

        t = np.zeros(n_sample)
        t[1:] = t_max * np.random.rand(n_sample - 1)

        y = np.sin((t - t_dead) * omega + phi)
        y += rel_noise * np.random.normal(size=n_sample)
        y *= amp
        y += offset

        # hold constant during dead time
        y = np.where(t > t_dead, y, offset + amp * np.sin(phi))

        # fix these fit parameters to a specific value
        const_dict = {
            't_dead': t_dead
            # 'phi': phi,
            # 'c': offset
        }
        p, p_err = sinusoid.fit(t,
                                y,
                                y_err=np.ones(y.shape) * amp * rel_noise,
                                evaluate_function=False,
                                evaluate_x_limit=[0, t_max],
                                constants=const_dict)
        print(p)
        print(p_err)

        print(p['omega'], "pm", 4 * p_err['omega'])
        print("fit =", p['a'], "pm", 4 * p_err['a'])
        print("fit =", p['c'], "pm", 4 * p_err['c'])
        print("fit =", p['phi'] % (2 * np.pi), "pm", 4 * p_err['phi'])
        print("fit lower bound =", p['t_dead'] - 4 * p_err['t_dead'])

        return p, p_err


    def test_random_data_fixed_phi(self):
        n_sample = 40
        t_max = 1e-6
        amp = 0.1
        rel_noise = 0.1
        offset = 40
        phi = np.pi
        t_dead = 0.0
        omega = np.pi / t_max * 9.3

        t = np.linspace(0, t_max, n_sample)

        y = np.sin((t - t_dead) * omega + phi)
        y += rel_noise * np.random.normal(size=n_sample)
        y *= amp
        y += offset

        # hold constant during dead time
        y = np.where(t > t_dead, y, offset + amp * np.sin(phi))

        # fix these fit parameters to a specific value
        const_dict = {
            # 't_dead': t_dead,
            'phi': phi,
            # 'c': offset
        }
        p, p_err = sinusoid.fit(t,
                                y,
                                y_err=np.ones(y.shape) * amp * rel_noise,
                                evaluate_function=False,
                                evaluate_x_limit=[0, t_max],
                                constants=const_dict)
        print(p)
        print(p_err)

        # self.assertAlmostEqual(omega, p['omega'], delta=4 * p_err['omega'])
        # self.assertAlmostEqual(amp, p['a'], delta=4 * p_err['a'])
        # self.assertAlmostEqual(offset, p['c'], delta=4 * p_err['c'])
        # self.assertAlmostEqual(phi % (2 * np.pi),
        #                        p['phi'] % (2 * np.pi),
        #                        delta=4 * p_err['phi'])
        # self.assertGreaterEqual(t_dead, p['t_dead'] - 4 * p_err['t_dead'])

    def test_pi_pulse(self):
        n_sample = 21
        t_max = 1e-4
        amp = 0.5
        rel_noise = 0.1
        offset = 0.5
        phi = np.pi / 2
        t_dead = 0.0
        omega = 1.3 * np.pi / t_max

        t = np.linspace(0, t_max, n_sample)

        y = np.sin((t - t_dead) * omega + phi)
        y += rel_noise * np.random.normal(size=n_sample)
        y *= amp
        y += offset

        # hold constant during dead time
        y = np.where(t > t_dead, y, offset + amp * np.sin(phi))

        # fix these fit parameters to a specific value
        const_dict = {
            't_dead': t_dead,
            # 'phi': phi,
            # 'c': offset
        }
        p, p_err = sinusoid.fit(t,
                                y,
                                y_err=np.ones(y.shape) * amp * rel_noise,
                                evaluate_function=False,
                                evaluate_x_limit=[0, t_max],
                                constants=const_dict)

        print(p)

    def test_pi_pulse_with_t_dead(self):
        n_sample = 21
        t_max = 1e-2
        amp = 0.5
        rel_noise = 0.05
        offset = 0.5
        phi = np.pi / 2
        t_dead = t_max * 0.24
        omega = 1.3 * np.pi / (t_max - t_dead)

        t = np.linspace(0, t_max, n_sample)

        y = np.sin((t - t_dead) * omega + phi)
        y += rel_noise * np.random.normal(size=n_sample)
        y *= amp
        y += offset

        # hold constant during dead time
        y = np.where(t > t_dead, y, offset + amp * np.sin(phi))

        # fix these fit parameters to a specific value
        const_dict = {
            # 't_dead': t_dead,
            'phi': phi,
            # 'c': offset
        }
        p, p_err = sinusoid.fit(t,
                                y,
                                y_err=np.ones(y.shape) * amp * rel_noise,
                                evaluate_function=False,
                                evaluate_x_limit=[0, t_max],
                                constants=const_dict)

        print(p)

    def test_delayed_split_data_with_t_dead(self):
        n0 = 40
        n1 = 30
        t0 = 1e-2
        t1 = t0 + 12.3e-4
        t2 = 2.1e-2
        t3 = t2 + 5.2e-4

        amp = 1.5
        rel_noise = 0.1
        offset = -5
        phi = np.pi / 2
        t_dead = (t3 - t2) * 0.14 + t0
        omega = np.pi / (t3 - t2) * 4.2

        t = np.concatenate((np.linspace(t0, t1, n0), np.linspace(t2, t3, n1)), axis=0)

        y = np.sin((t - t_dead) * omega + phi)
        y += rel_noise * np.random.normal(size=n0 + n1)
        y *= amp
        y += offset

        # hold constant during dead time
        y = np.where(t > t_dead, y, offset + amp * np.sin(phi))

        # fix these fit parameters to a specific value
        const_dict = {
            # 't_dead': t_dead,
            'phi': phi,
            # 'c': offset
        }
        p, p_err = sinusoid.fit(t,
                                y,
                                y_err=np.ones(y.shape) * amp * rel_noise,
                                evaluate_function=False,
                                evaluate_x_limit=[0, t3],
                                constants=const_dict)

        print(t)
        print(p)












#
# class SinusoidTest(EnvExperiment):
#     def build(self):
#         self.setattr_device("core")  # not strictly needed here
#         # optional: seed for reproducibility
#         self.setattr_argument("seed", NumberValue(0, step=1))
#         # simple logging
#         logging.basicConfig(level=logging.INFO)
#         self.logger = logging.getLogger(__name__)
#
#     def prepare(self):
#         np.random.seed(int(self.seed))
#
#         # ---- parameters (from your test_random_data) ----
#         self.n_sample = 30
#         self.t_max    = 1e-1
#         self.amp      = 1.5
#         self.rel_noise= 0.1
#         self.offset   = 40.0
#         self.phi      = np.pi/2
#         self.t_dead   = 0.0
#         self.omega    = np.pi/self.t_max * 8.2
#
#         # ---- generate synthetic data on host ----
#         t = np.zeros(self.n_sample)
#         t[1:] = self.t_max * np.random.rand(self.n_sample - 1)
#
#         y = np.sin((t - self.t_dead) * self.omega + self.phi)
#         y += self.rel_noise * np.random.normal(size=self.n_sample)
#         y *= self.amp
#         y += self.offset
#         # hold constant during dead time
#         y = np.where(t > self.t_dead, y, self.offset + self.amp * np.sin(self.phi))
#
#         self.t = t
#         self.y = y
#
#     def run(self):
#         # Do the fit on the host
#         const_dict = {"t_dead": self.t_dead}
#         p, p_err = sinusoid.fit(
#             self.t,
#             self.y,
#             y_err=np.ones_like(self.y) * self.amp * self.rel_noise,
#             evaluate_function=False,
#             evaluate_x_limit=[0, self.t_max],
#             constants=const_dict
#         )
#
#         # ---- Publish everything to datasets (visible in Dashboard) ----
#         # raw data
#         self.set_dataset("fit_data.t", self.t, broadcast=True)
#         self.set_dataset("fit_data.y", self.y, broadcast=True)
#
#         # ground truth (known) params
#         self.set_dataset("fit_truth.omega",  self.omega,  broadcast=True)
#         self.set_dataset("fit_truth.a",      self.amp,    broadcast=True)
#         self.set_dataset("fit_truth.c",      self.offset, broadcast=True)
#         self.set_dataset("fit_truth.phi",    float(np.mod(self.phi, 2*np.pi)), broadcast=True)
#         self.set_dataset("fit_truth.t_dead", self.t_dead, broadcast=True)
#
#         # fitted params + 4σ errors
#         for k, v in p.items():
#             self.set_dataset(f"fit.{k}", float(v), broadcast=True)
#         for k, v in p_err.items():
#             self.set_dataset(f"fit_err4.{k}", float(4*v), broadcast=True)
#
#         # convenience: modulo 2π for phase
#         phi_fit_mod  = float(np.mod(p["phi"], 2*np.pi))
#         phi_true_mod = float(np.mod(self.phi, 2*np.pi))
#         self.set_dataset("fit.phi_mod_2pi",  phi_fit_mod,  broadcast=True)
#         self.set_dataset("fit_truth.phi_mod_2pi", phi_true_mod, broadcast=True)
#
#         # Optional: also log to controller log (may show in Dashboard log panel)
#         self.logger.info("omega: true=%.6g fit=%.6g ± %.6g",
#                          self.omega, p["omega"], 4*p_err["omega"])
#         self.logger.info("a:     true=%.6g fit=%.6g ± %.6g",
#                          self.amp, p["a"], 4*p_err["a"])
#         self.logger.info("c:     true=%.6g fit=%.6g ± %.6g",
#                          self.offset, p["c"], 4*p_err["c"])
#         self.logger.info("phi:   true=%.6g fit=%.6g ± %.6g",
#                          phi_true_mod, phi_fit_mod, 4*p_err["phi"])
#         self.logger.info("t_dead: true=%.6g fit lower bound=%.6g",
#                          self.t_dead, p["t_dead"] - 4*p_err["t_dead"])
#
#         # Add a quick pass/fail flag like a test would
#         ok_omega = abs(p["omega"] - self.omega) <= 4*p_err["omega"]
#         ok_a     = abs(p["a"]     - self.amp)   <= 4*p_err["a"]
#         ok_c     = abs(p["c"]     - self.offset)<= 4*p_err["c"]
#         ok_phi   = abs(phi_fit_mod - phi_true_mod) <= 4*p_err["phi"]
#         ok_tdead = (self.t_dead >= p["t_dead"] - 4*p_err["t_dead"])
#
#         self.set_dataset("fit_check.ok_omega", ok_omega, broadcast=True)
#         self.set_dataset("fit_check.ok_a",     ok_a,     broadcast=True)
#         self.set_dataset("fit_check.ok_c",     ok_c,     broadcast=True)
#         self.set_dataset("fit_check.ok_phi",   ok_phi,   broadcast=True)
#         self.set_dataset("fit_check.ok_tdead", ok_tdead, broadcast=True)
#
#         # Optional: single summary string
#         summary = (
#             f"omega true={self.omega:.6g} fit={p['omega']:.6g} ±{4*p_err['omega']:.2g} | "
#             f"a true={self.amp:.6g} fit={p['a']:.6g} ±{4*p_err['a']:.2g} | "
#             f"c true={self.offset:.6g} fit={p['c']:.6g} ±{4*p_err['c']:.2g} | "
#             f"phi true(mod 2π)={np.mod(self.phi,2*np.pi):.6g} fit(mod 2π)={phi_fit_mod:.6g} ±{4*p_err['phi']:.2g} | "
#             f"t_dead check: true={self.t_dead:.6g} >= {p['t_dead'] - 4*p_err['t_dead']:.6g}"
#         )
#
#         self.set_dataset("fit_summary", summary, broadcast=True)