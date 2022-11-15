from utils import *
from tqdm import tqdm
from scipy.optimize import newton

U_guess = -.311*eV                                  # correcto -0.311 *eV
E_ = np.linspace(E_C, 50*eV, 50000)
T = 300
V_DS = 0.3
V_G = 0.4
args = (T, V_DS, V_G)

U_ana = deme_un_U_mi_Rey(V_G_array)

# Punto fijo
# U_array = np.linspace(-0.5*eV, -0.2*eV, 2000)
# U_fixed = np.array([get_U_fixed_point(U_array, T, V_DS, VG) for VG in tqdm(V_G_array)])
# plt.plot(V_G_array, U_fixed/eV, label='Fixed point')
# plt.plot(V_G_array, U_ana/eV, label='Analytical')
# plt.show()

# Root finding
# U_array = np.linspace(-0.5*eV, -0.1*eV, 100)
# def root_alg(U, *args):
#     return (iter_alg(U, *args)[2] - U) / eV

# U_new = np.array([newton(root_alg, U, args=args) for U in U_array])

# plt.plot(U_array/eV, U_new/eV)
# plt.plot(U_array/eV, U_array/eV)
# plt.show()


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

print(f'N0 = {N0:.3e}, N = {N:.3e}')
print(f"U = {U_new / eV:.6f} eV")

fig, ax = plt.subplots(1, 3, figsize=(15, 5))
ax[0].plot(N_array, '.')
ax[0].set_title("N")
ax[1].plot(N0_array, '.')
ax[1].set_title("N0")
ax[2].plot(U_array, '.')
ax[2].set_title("U")
plt.show()