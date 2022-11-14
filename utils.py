import numpy as np
from scipy.constants import k, e, hbar, eV
import matplotlib.pyplot as plt

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
mu_s = E_F
mu_d = lambda _V_DS: E_F - e*_V_DS          # J

C_Q = e**2 *m_eff*A/(2*np.pi*hbar**2)


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


def integrate(integrand, x):
    """
    Esta función integra cualquier función capaz de recibir un escalar o np.array y la integra en los
    límites especificados usando el método de trapecios
    :param integrand: es el objeto de python que representa la función a integrar
    :param min_value: es el límite inf de la integral
    :param max_value: es el lim sup de la integral
    :return: devuelve un escalar que es el valor de la integral
    """
    return np.trapz(integrand, x=x)

def calculate_N(E: np.ndarray|float, T, U, V_DS):
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
    integrand1 = g(E-E_C-U) * f(E - mu_s, T)            # N_s
    integrand2 = g(E-E_C-U) * f(E - mu_d(V_DS), T)      # N_d
    integrand = (integrand1 + integrand2) / 2
    #plt.plot(E/eV, integrand1, label='N_s')
    #plt.plot(E/eV, integrand2, label='N_d')
    #plt.plot(E/eV, integrand, label='N')
    #plt.vlines(E_C/eV, 0, 1, colors='r', linestyles='dashed', label='E_C')
    #plt.vlines(E_C/eV+U/eV, 0, 1, colors='g', linestyles='dashed', label='E_C+U')
    #plt.legend()
    #plt.show()
    #breakpoint()
    return integrate(integrand, x=E) * A * m_eff / (np.pi * hbar**2)


def calculate_N0(E: np.ndarray|float, T):
    return integrate(g(E - E_C) * f(E - E_F, T), x=E) * A * m_eff / (np.pi * hbar**2)


def calculate_U(N, N0, V_G):
    """
    Esta función calcula el valor de U usando el valor de N calculado anteriormente
    la idea es poder comparar el U que salga de acá con la guess
    :param N: Número de electrones
    :return: Un potencial U actualizado para ser más consistente con N
    """
    U_es = -V_G * e
    U_C = e**2 * (N - N0) / C_es            # Es muy positivo :(
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

def iter_alg(U, T, V_DS, V_G):
    """
    Esta es una iteración del algoritmo, recibe una guess y devuelve un U nuevo
    :param U: U en la iteración anterior
    :return: U nuevo
    """
    Emin = (E_C + U)
    Emax = -2*eV
    E_N = np.linspace(Emin, Emax, 10000)
    E_N0 = np.linspace(E_C, Emax, 10000)
    N = calculate_N(E_N, T, U, V_DS)
    N0 = calculate_N0(E_N0, T)
    U_new = calculate_U(N, N0, V_G)
    return N, N0, U_new

def calculate_I(E: np.ndarray|float, U, T, V_DS):
    """Funciona, no cambiar, gracias."""
    integrand1 = np.sqrt(2*m_eff*(E-E_C-U), where=E-E_C-U>=0) * (E-E_C-U >= 0)
    integrand2 =  f(E - mu_s, T) - f(E - mu_d(V_DS), T)
    integrand = integrand1 * integrand2
    return integrate(integrand, x=E) * e * Width / (np.pi * hbar)**2

def deme_un_U_mi_Rey(V_G):
    """Ojooooooo válido pa 0 K"""
    eta_0=1
    V_T=(E_C-mu_s)/(eta_0*e)
    C_Q=e**2*m_eff*Width*Length/(2*np.pi*hbar**2)
    eta=C_G/(C_G+C_Q)
    U=-eta*e*(V_G-V_T)-eta_0*e*V_T
    return U
