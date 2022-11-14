from utils import *
from tqdm import tqdm

U_guess = -0.311*eV                                  # correcto -0.311 *eV
E_ = np.linspace(E_C, 50*eV, 50000)
T = 298
V_DS = 0.
V_G = .4
args = (T, V_DS, V_G)

U_array = deme_un_U_mi_Rey(V_G_array)

#for i, VG in enumerate(V_G_array):
#    plt.plot(V_DS_array, [calculate_I(E_+U_array[i], U_array[i], T, VDS) for VDS in V_DS_array])
#plt.plot(V_DS_array, [calculate_I(E_, U_guess, T, VDS) for VDS in V_DS_array])
#plt.show()
#print(calculate_I(E_, U_guess, T, V_DS))

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
ax[1].plot(N0_array, '.')
ax[1].set_title("N0")
ax[2].plot(U_array, '.')
ax[2].set_title("U")
plt.show()
