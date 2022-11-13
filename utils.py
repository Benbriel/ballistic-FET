import numpy as np
from scipy.integrate import quad

# constants
mu_s = 1
mu_d = 1
N0 = 0
q = 1
V_gs = 1
C_g = 1
E_min = 1
E_max = 1
k_B = 8.617333262 * 1e-5    # ev/K
T = 1                       # K
def g(E):
    """
    Esta función calcula la densidad de estados para una energía en particular
    :param E: Es la energía para la cual se quiere calcular la función de densidad de estados
    :return: devuelve un np.array con los resultados de la entrada que se le da
    """
    pass


def f(E, mu):
    """
    Esta es la distribución Fermi-Dirac
    :param E: Energía en eV
    :param mu: Potencial químico en eV
    :return: valor de la función FD np.array
    """
    return 1 / (np.exp((E - mu)/(k_B * T)) + 1)



def integrate(integrand, min_value, max_value):
    """
    Esta función integra cualquier función capaz de recibir un escalar o np.array y la integra en los
    límites especificados
    :param integrand: es el objeto de python que representa la función a integrar
    :param min_value: es el límite inf de la integral
    :param max_value: es el lim sup de la integral
    :return: devuelve un escalar que es el valor de la integral
    """
    return quad(integrand, min_value, max_value)

def calculate_N(E_min, E_max, U, mu_s=mu_s, mu_d=mu_d):
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
    integrand1 = lambda E: g(E-U) * f(E, mu_s)
    integrand2 = lambda E: g(E-U) * f(E, mu_d)
    integral1 = integrate(integrand1, E_min, E_max)
    integral2 = integrate(integrand2, E_min, E_max)
    return (integral1 + integral2) / 2


def calculate_N0():
    pass

def calculate_U(N):
    """
    Esta función calcula el valor de U usando el valor de N calculado anteriormente
    la idea es poder comparar el U que salga de acá con la guess
    :param N: Número de electrones
    :return: Un potencial U actualizado para ser más consistente con N
    """
    return -q * V_gs + q ** 2 * (N - N0) / C_g


def convergence_criteria(value_before, value_now, th):
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

def iter_alg(U):
    """
    Esta es una iteración del algoritmo, recibe una guess y devuelve un U nuevo
    :param U: U en la iteración anterior
    :return: U nuevo
    """
    N = calculate_N(E_min, E_max, U)
    U_new = calculate_U(N)
    return U_new
