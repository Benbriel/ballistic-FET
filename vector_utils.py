import numpy as np
from scipy.constants import k, e, hbar, eV
import matplotlib.pyplot as plt
from tqdm import tqdm
import os

# Define Constants
nm = 1e-9                                   # m
Length = 40 * nm                            # m
Width = 3 * Length                          # m
A = Length * Width                          # m^2
E_F = -5. * eV                              # J
E_C = -4.7 * eV                             # J
# tau_d = tau_s -> 0. tal que tau_d / tau_s -> 1
C_es = C_G = 0.1 * 1e-15                    # F
T_array = np.array([1., 298])               # K
V_DS_array = np.linspace(0, .5, 100)        # V
V_G_array = np.array([.3, .35, .4, .45, .5])# V
m_eff = 9.1 / 2 * 1e-31                     # kg

# Algunas variables
beta = lambda T: 1/(k*T)                    # J^-1
mu_s = E_F                                  # J
mu_d = lambda _V_DS: E_F - e*_V_DS          # J
alpha = A * m_eff / (np.pi * hbar**2)       # J^-1
C_Q = e**2 * alpha / 2
eta = C_G / (C_G + C_Q)
eta_0 = 1
V_T = (E_C - mu_s) / (eta_0 * e)
E_ = lambda Emin: np.linspace(Emin, 0*eV, 1000)

def g(E):
    """
    Esta función calcula la densidad de estados para una energía en particular
    :param E: Es la energía para la cual se quiere calcular la función de densidad de estados
    :return: devuelve un np.array con los resultados de la entrada que se le da
    """
    return (E >= 0)

def f(E: np.ndarray, T: np.ndarray):
    """
    Esta es la distribución Fermi-Dirac
    :param E: Energía en eV
    :param mu: Potencial químico en eV
    :return: valor de la función FD np.array
    """
    expo = beta(T) * E
    ans = np.empty_like(expo)
    is_neg = expo < 0; is_pos = ~is_neg
    ans[is_neg] = 1 / (np.exp(expo[is_neg]) + 1)
    ans[is_pos] = np.exp(-expo[is_pos]) / (np.exp(-expo[is_pos]) + 1)
    return ans

def get_N(U: np.ndarray, T: np.ndarray, V_DS: np.ndarray) -> np.ndarray:
    """
    Esta función calcula el valor de N según la ecuación (5.13), (5.15) y (5.16) del libro del curso
    :param E_min: la energía mínima por donde partir la integral, en rigor es menos infinito pero
    por temas numéricos es recomendable solo usar valores para los cuales el integrando empieza
    a ser despreciable
    :param E_max: la energía mínima por donde partir la integral, en rigor es infinito pero
    por temas numéricos es recomendable solo usar valores para los cuales el integrando empieza
    :param U: El calor que en el algoritmo sale como 'U guess' para el potencial U
    :param mu_s: potencial químico source en eV
    :param mu_d: potencial químido drain en eV
    :return: devuelve el valor de N calculado a partir de estadisticas Fermi-Dirac
    """
    # integrand.shape == (len(E), len(T), len(V_DS), len(V_G))
    E = E_(E_C+U)
    T = T[None, :, None, None]
    V_DS = V_DS[None, None, :, None]
    integrand1 = (g(E-E_C-U) * f(E - mu_s, T))
    integrand2 = g(E-E_C-U) * f(E - mu_d(V_DS), T)
    integrand = (integrand1 + integrand2) / 2
    return np.trapz(integrand, x=E, axis=0) * alpha

def get_N0(T):
    E = E_(E_C)[:, None]
    T = T[None, :]
    integrand = g(E-E_C) * f(E - mu_s, T)
    return np.trapz(integrand, x=E, axis=0)[:, None, None] * alpha

def get_U(N, N0, V_G):
    """
    Esta función calcula el valor de U usando el valor de N calculado anteriormente
    la idea es poder comparar el U que salga de acá con la guess
    :param N: Número de electrones
    :return: Un potencial U actualizado para ser más consistente con N
    """
    U_es = -V_G * e
    U_C = e**2 * (N - N0) / C_es
    return  (U_es + U_C)

def get_new_U(U_guess, T, V_DS, V_G, delta=0.005):
    """
    Esta es una iteración del algoritmo, recibe una guess y devuelve un U nuevo
    :param U: U en la iteración anterior
    :return: U nuevo
    """
    N = get_N(U_guess, T, V_DS)
    N0 = get_N0(T)
    U = get_U(N, N0, V_G)
    return U_guess + delta * (U - U_guess)

def iter_U(U_guess: float|np.ndarray, T: np.ndarray,
           V_DS: np.ndarray, V_G: np.ndarray, n_iter=2000, delta=0.005):
    """
    Esta función itera el algoritmo para obtener un U consistente con N
    :param U_guess: U inicial
    :param n_iter: número de iteraciones
    :return: U final
    """
    U_iter = np.empty((n_iter, T.size, V_DS.size, V_G.size))
    if isinstance(U_guess, float):
        U = np.array([U_guess])[:, None, None, None]
    else:
        assert U_guess.shape == (1, T.size, V_DS.size, V_G.size)
        U = U_guess
    for i in tqdm(range(n_iter)):
        U = get_new_U(U, T, V_DS, V_G, delta=delta)
        U_iter[i] = U[0]
    return U_iter

def get_I(U, T, V_DS):
    """
    Esta función calcula la corriente I
    :param U: potencial U
    :param T: temperatura
    :param V_DS: voltaje drain-source
    :return: corriente I
    """
    E = E_(E_C+U)
    T = T[None, :, None, None]
    V_DS = V_DS[None, None, :, None]
    integrand1 = np.sqrt(2*m_eff*(E-E_C-U), where=E-E_C-U>=0) * g(E-E_C-U)
    integrand2 = f(E - mu_s, T) - f(E - mu_d(V_DS), T)
    integrand = integrand1 * integrand2
    I = np.trapz(integrand, x=E, axis=0) * e * Width / (np.pi * hbar)**2
    return I

def get_I_0K(V_DS: np.ndarray, V_G: np.ndarray):
    V_DS = V_DS[None, None, :, None]
    V_G = V_G[None, None, None, :]
    const = (e*Width)/(np.pi**2*hbar**2)*np.sqrt(8*m_eff/9)*((eta*e)**(3/2))
    linear = (V_G-V_T)**(3/2) - (V_G-V_T-V_DS/eta) * np.sqrt(V_G-V_T-V_DS/eta, where=V_G-V_T-V_DS/eta>=0)
    sat = (V_G-V_T)**(3/2)
    is_linear = (V_DS <= eta*(V_G - V_T))
    return const * (is_linear * linear + (~is_linear) * sat)

def plot_U(U, U_iter):
    for vg in range(len(V_G_array)):
        plt.plot(V_DS_array, U[0, 0, :, vg]/eV, label=f'$V_G$ = {V_G_array[vg]} V')
    plt.legend()
    plt.xlabel('V_DS [V]')
    plt.ylabel('U [eV]')
    plt.tight_layout()
    plt.show()

    plt.plot(U_iter[:, 0, -1, -2]/eV, label='T = 1 K')
    plt.plot(U_iter[:, 1, -1, -2]/eV, label='T = 298 K')
    plt.legend()
    plt.xlabel('iteración')
    plt.ylabel('U [eV]')
    plt.tight_layout()
    plt.show()
