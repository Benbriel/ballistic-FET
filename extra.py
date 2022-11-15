from utils import *
from tqdm import tqdm

# Constantes
U_guess = (-.29+.5)*eV                                  # correcto -0.311 *eV
E_ = np.linspace(E_C, 50*eV, 50000)
T = 300
V_DS = 0.3
V_G = 0.4
args = (T, V_DS, V_G)

U_ana = deme_un_U_mi_Rey(V_G_array)

#for i, VG in enumerate(V_G_array):
#    plt.plot(V_DS_array, [calculate_I(E_+U_array[i], U_array[i], T, VDS) for VDS in V_DS_array])
#plt.plot(V_DS_array, [calculate_I(E_, U_guess, T, VDS) for VDS in V_DS_array])
#plt.show()
#print(calculate_I(E_, U_guess, T, V_DS))

# Punto fijo
U_array = np.linspace(-0.5*eV, -0.2*eV, 2000)
U_fixed = np.array([get_U_fixed_point(U_array, T, V_DS, VG) for VG in tqdm(V_G_array)])
plt.plot(V_G_array, U_fixed/eV, label='Fixed point')
plt.plot(V_G_array, U_ana/eV, label='Analytical')
plt.show()

# Puntos fijos
U_aux_array = np.linspace(-0.5*eV, 0.*eV, 1000)
Nnew, N0new, Unew = np.array([iter_alg(U, T, V_DS, V_G) for U in U_aux_array]).T
Nnew2, N0new2, Unew2 = np.array([iter_alg(U, T, V_DS, V_G) for U in Unew]).T
Nnew3, N0new3, Unew3 = np.array([iter_alg(U, T, V_DS, V_G) for U in Unew2]).T
Nnew4, N0new4, Unew4 = np.array([iter_alg(U, T, V_DS, V_G) for U in Unew3]).T
plt.plot(U_aux_array/eV, U_aux_array/eV)
plt.plot(U_aux_array/eV, Unew/eV)
plt.plot(U_aux_array/eV, Unew2/eV, c='r')
plt.plot(U_aux_array/eV, Unew3/eV, c='g')
plt.plot(U_aux_array/eV, Unew4/eV, c='y')
plt.show()


# Iteraciones
n_iter = 1000
N_array = np.zeros(n_iter)
N0_array = np.zeros(n_iter)
U_array = np.zeros(n_iter + 1)
U_array[0] = U_guess / eV
for i in tqdm(range(n_iter)):
    N, N0, U_new = iter_alg(U_guess, *args)
    U_guess = U_new
    N_array[i] = N
    N0_array[i] = N0
    U_array[i+1] = U_new / eV
    #breakpoint()

print(f'N0 = {N0:.3e}, N = {N:.3e}')
print(f"U = {U_new / eV:.6f} eV")

fig, ax = plt.subplots(1, 3, figsize=(15, 5))
ax[0].plot(N_array, '.')
ax[0].set_title("N")
#ax[0].set_yscale('log')
ax[1].plot(N0_array, '.')
ax[1].set_title("N0")
#ax[1].set_yscale('log')
ax[2].plot(U_array, '.')
ax[2].set_title("U")
#ax[2].set_yscale('log')
plt.show()