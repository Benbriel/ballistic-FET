from scipy.constants import k, e, hbar, eV
from scipy.integrate import quad
import matplotlib.pyplot as plt
import numpy as np
from utils import *

# Define Constants
hbar /= eV
e /= eV
A = 1.
E_F = -5. * eV
mu_s = 1
mu_d = 1
tau_d = 1.
tau_s = 1.
# T = 1
N0 = 1
C_es = 1
A = 1.
m = 1.
infty = 100.
beta = lambda T: 1/(k*T) * eV

def f(E, mu, T):
    return 1 / (1 + np.exp(np.clip((E - mu) * beta(T), -100, 1e2)))


def U_potential(N, q=1):
    return q**2 *  (N - N0) / C_es


def g(E, E_n, nmax=0):
    return A*m/(np.pi*hbar)*(E > E_n)


def calc_N(U, T):
    e_dist = lambda E: (tau_d * f(E, mu_s, T) + tau_s * f(E, mu_d, T)) / (tau_s + tau_d)
    g_dist = lambda E: g(E - U, 0) * e_dist(E)
    integrand = lambda E: g_dist(E)
    E_plot = np.linspace(-20, 20, 1000)
    plt.plot(E_plot, integrand(E_plot))
    plt.show()
    breakpoint()
    integral = quad(integrand, -infty, infty)[0]
    return integral


def get_U_optimized(U_guess, T, q, th=1e-4, max_iter=10):
    if max_iter == 0:
        breakpoint()
        raise ValueError("Max Iterations Reached")
    N = calc_N(U_guess, T)
    U_new = U_potential(N)
    if np.abs(U_new - U_guess) < th:
        return U_new
    return get_U_optimized(U_new, T, q, th, max_iter-1)


def calc_I(U, T, q):
    e_dist = lambda E: (f(E, mu_s, T) - f(E, mu_d, T)) / (tau_s + tau_d)
    g_dist = lambda E: g(E - U, 0)
    integral = quad(lambda E: e_dist(E) * g_dist(E), -infty, infty)[0]
    return integral


if __name__ == "__main__":
    T = 1
    U = get_U_optimized(1., T, e)
    print(U)
    print(calc_I(U, T, e))