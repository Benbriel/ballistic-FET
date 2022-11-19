import numpy as np
from scipy.constants import k, e, hbar, eV
import matplotlib.pyplot as plt
from tqdm import tqdm

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
E_ = lambda Emin: np.linspace(Emin, 0*eV, 10000)


def g(E):
    """
    Esta función calcula la densidad de estados para una energía en particular
    :param E: Es la energía para la cual se quiere calcular la función de densidad de estados
    :return: devuelve un np.array con los resultados de la entrada que se le da
    """
    return (E >= 0)


def f(E, T):
    """
    Esta es la distribución Fermi-Dirac
    :param E: Energía en eV
    :param mu: Potencial químico en eV
    :return: valor de la función FD np.array
    """
    return np.piecewise(E, [E < 0, E >= 0],
    [
        lambda E: 1 / (np.exp(beta(T) * E) + 1),
        lambda E: np.exp(-beta(T) * E) / (np.exp(-beta(T) * E) + 1)
        ]
    )

def calculate_N(T, U, V_DS: float):
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
    E = E_(E_C+U)
    integrand1 = g(E-E_C-U) * f(E - mu_s, T)            # N_s
    integrand2 = g(E-E_C-U) * f(E - mu_d(V_DS), T)      # N_d
    integrand = (integrand1 + integrand2) / 2
    return np.trapz(integrand, x=E) * alpha

def N_ana(U, T, V_DS):
    return alpha * (
        (np.log(np.exp(beta(T)*mu_s) + np.exp(beta(T)*(E_C+U))) + \
            np.log(np.exp(beta(T)*mu_d(V_DS)) + np.exp(beta(T)*(E_C+U)))) / \
                (2*beta(T)) - (E_C+U)
    )

def N0_ana(T):
    return alpha * (np.log(np.exp(beta(T)*E_F) + np.exp(beta(T)*E_C)) / beta(T) - E_C)


def calculate_N0(T):
    E = E_(E_C)
    return np.trapz(g(E - E_C) * f(E - E_F, T), x=E) * alpha


def calculate_U(N, N0, V_G):
    """
    Esta función calcula el valor de U usando el valor de N calculado anteriormente
    la idea es poder comparar el U que salga de acá con la guess
    :param N: Número de electrones
    :return: Un potencial U actualizado para ser más consistente con N
    """
    U_es = -V_G * e
    U_C = e**2 * (N - N0) / C_es
    return  (U_es + U_C)


def convergence_criteria(value_before, value_now, th=1e-6):
    """
    Este es el criterio de convergencia usado para evaluar si el main loop debe detenerse
    :param value_before: es el valor de interés en la iteración anterior
    :param value_now: es el valor de interés en la iteración actual
    :param th: es el umbral de la diferencia para evaluar convergencia
    :return: devuelve True si hay convergencia y False en caso contrario
    """
    return np.abs(value_before - value_now) < th

def iter_alg(U, T, V_DS, V_G, delta=0.005):
    """
    Esta es una iteración del algoritmo, recibe una guess y devuelve un U nuevo
    :param U: U en la iteración anterior
    :return: U nuevo
    """
    N = calculate_N(T, U, V_DS)
    N0 = calculate_N0(T)
    U_new = calculate_U(N, N0, V_G)
    return U + delta * (U_new - U)

def calculate_I(U, T, V_DS):
    """Funciona, no cambiar, gracias."""
    E = E_(E_C+U)
    integrand1 = np.sqrt(2*m_eff*(E-E_C-U), where=E-E_C-U>=0) * (E-E_C-U >= 0)
    integrand2 =  f(E - mu_s, T) - f(E - mu_d(V_DS), T)
    integrand = integrand1 * integrand2
    return np.trapz(integrand, x=E) * e * Width / (np.pi * hbar)**2

def deme_un_U_mi_Rey(V_G):
    """Ojooooooo válido pa 0 K"""
    eta_0=1
    V_T=(E_C-mu_s)/(eta_0*e)
    C_Q=e**2*m_eff*Width*Length/(2*np.pi*hbar**2)
    eta=C_G/(C_G+C_Q)
    U=-eta*e*(V_G-V_T)-eta_0*e*V_T
    return U

def get_U_fixed_point(U_array: np.ndarray, T, V_DS, V_G):
    Unew = np.array([iter_alg(U, T, V_DS, V_G) for U in U_array]).T
    is_zero = np.diff(np.sign(Unew - U_array)) != 0
    return Unew[1:][is_zero][0]

def get_U_iterative(U, T, V_DS, V_G, niter=2000) -> np.ndarray:
    U_array = np.zeros(niter)
    U_array[0] = U
    for i in range(1, niter):
        U = iter_alg(U, T, V_DS, V_G)
        U_array[i] = U

    return U_array
