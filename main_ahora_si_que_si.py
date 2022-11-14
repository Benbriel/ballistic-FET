from scipy.constants import k, e, hbar, eV
from scipy.integrate import quad
import matplotlib.pyplot as plt
import numpy as np

# Define Constants
hbar /= eV                                  # eV*s
nm = 1e-9                                   # m
Length = 40 * nm                            # m
Width = 3 * Length                          # m
A = Length * Width                          # m^2
E_F = -5.                                   # eV
E_C = -4.7                                  # eV
# tau_d = tau_s -> 0. tal que tau_d / tau_s -> 1
C_es = C_G = 0.1 * 1e-15                    # F
T_array = np.array([1., 298])               # K
V_DS_array = np.linspace(0, .5, 100)        # V
V_G_array = np.array([.3, .35, .4, .45, .5])# V
m_eff = 9.1 / 2 * 1e-31                     # kg

# Algunas variables
beta = lambda T: 1/(k*T) * eV               # eV^-1
mu_s = E_F
mu_d = lambda _V_DS: E_F + e*_V_DS / eV     # eV

# Para testear
E_plot = np.linspace(-6, -2, 100000, dtype=np.longdouble)

def f(E, T):
    return 1 / (1 + np.exp(np.clip(E * beta(T), -np.inf, 1e2)))

def f(E, T):
    neg_exp = np.exp(-E * beta(T))
    return neg_exp / (1 + neg_exp)

def g(E):
    E *= (E > 0)                            # step function
    return A / (np.pi * hbar)**2

def get_N(T, mu, E_n=E_C):
    integrand = lambda E: g(E - E_n) * f(E - mu, T)
    integral = quad(integrand, E_n, 500*np.abs(E_C))[0]
    return integral

def U_potential(N, q=e):
    return
    return q**2 *  (N - N0) / C_es / eV

def get_best_U(U_guess, T, q, th=1e-4, max_iter=10):
    return
    if max_iter == 0:
        breakpoint()
        raise ValueError("Max Iterations Reached")
    N = calc_N(U_guess, T, q)
    U_new = U_potential(N)
    if np.abs(U_new - U_guess) < th:
        return U_new
    return get_best_U(U_new, T, q, th, max_iter-1)

if __name__ == "__main__":
    plt.plot(E_plot, f(E_plot - mu_s, 150), label="f")
    plt.plot(E_plot, g(E_plot-E_C), label="g")
    plt.twinx()
    plt.plot(E_plot, f(E_plot - mu_s, 150)*g(E_plot-E_C), label="f*g", c='r')
    plt.legend()
    plt.show()

breakpoint()
