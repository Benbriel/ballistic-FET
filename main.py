from utils import *

U_guess = 0
T = 300
V_DS = 0
V_G = 0
# Emin = E_C + U_guess
Emax = 100 * np.abs(E_C)
args = (Emax, T, V_DS, V_G)

U_array = []
for i in range(20):
    U_new = iter_alg(U_guess, E_C + U_guess, *args)
    print('U_new =', U_new)
    U_guess = U_new
    U_array.append(U_new)

plt.plot(U_array)
plt.yscale('symlog')
plt.show()

breakpoint()

calculate_N(E_C+U_guess, Emax, 300, U_guess, V_DS_array[-1])

breakpoint()
