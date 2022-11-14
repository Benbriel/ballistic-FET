import numpy as np
from scipy.integrate import quad
from scipy.constants import k, e, hbar, eV
import matplotlib.pyplot as plt

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

def g(E):
    """
    Esta función calcula la densidad de estados para una energía en particular
    :param E: Es la energía para la cual se quiere calcular la función de densidad de estados
    :return: devuelve un np.array con los resultados de la entrada que se le da
    """
    return A * m_eff / (np.pi * hbar**2) * (E > 0)


def f(E, T):
    """
    Esta es la distribución Fermi-Dirac
    :param E: Energía en eV
    :param mu: Potencial químico en eV
    :return: valor de la función FD np.array
    """
    neg_exp = np.exp(-E * beta(T))
    return neg_exp / (1 + neg_exp)


def integrate(integrand, min_value, max_value):
    """
    Esta función integra cualquier función capaz de recibir un escalar o np.array y la integra en los
    límites especificados
    :param integrand: es el objeto de python que representa la función a integrar
    :param min_value: es el límite inf de la integral
    :param max_value: es el lim sup de la integral
    :return: devuelve un escalar que es el valor de la integral
    """
    return quad(integrand, min_value, max_value)[0]

def calculate_N(E_min, E_max, T, U, V_DS):
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
    integrand1 = lambda E: g(E-E_C-U) * f(E - mu_s, T)
    integrand2 = lambda E: g(E-E_C-U) * f(E - mu_d(V_DS), T)
    integral1 = integrate(integrand1, E_min, E_max)
    integral2 = integrate(integrand2, E_min, E_max)
    return (integral1 + integral2) / 2


def calculate_N0(T):
    return integrate(lambda E: g(E-E_C) * f(E - E_F, T), E_C, np.inf)


def calculate_U(N, N0, V_G):
    """
    Esta función calcula el valor de U usando el valor de N calculado anteriormente
    la idea es poder comparar el U que salga de acá con la guess
    :param N: Número de electrones
    :return: Un potencial U actualizado para ser más consistente con N
    """
    return -V_G + e ** 2 * (N - N0) / (C_es * eV)


def convergence_criteria(value_before, value_now, th=1e-6):
    """
    Este es el criterio de convergencia usado para evaluar si el main loop debe detenerse
    :param value_before: es el valor de interés en la iteración anterior
    :param value_now: es el valor de interés en la iteración actual
    :param th: es el umbral de la diferencia para evaluar convergencia
    :return: devuelve True si hay convergencia y False en caso contrario
    """
    if np.abs(value_before - value_now) < th:
        return True
    else:
        return False

def iter_alg(U, Emin, Emax, T, V_DS, V_G):
    """
    Esta es una iteración del algoritmo, recibe una guess y devuelve un U nuevo
    :param U: U en la iteración anterior
    :return: U nuevo
    """
    N = calculate_N(Emin, Emax, T, U, V_DS)
    print('N: ', N)
    N0 = calculate_N0(T)
    print('N0: ', N0)
    U_new = calculate_U(N, N0, V_G)
    return U_new
