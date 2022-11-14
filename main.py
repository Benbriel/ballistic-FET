from utils import *
from tqdm import tqdm

U_guess = 0 * eV
E_ = np.linspace(E_C + U_guess, 50*eV, 50000)
T = 300
V_DS = .5
V_G = .5
args = (T, V_DS, V_G)

U_array = np.linspace(-100, -1, 100) * eV
for VG in V_G_array:
    plt.plot(V_DS_array, [calculate_I(E_, -VG, T, -V_DS) for V_DS in V_DS_array])
plt.show()
#print(calculate_I(E_, U_guess, T, V_DS))
#breakpoint()
N_array = []
N0_array = []
U_array = []
for i in tqdm(range(1000)):
    N, N0, U_new = iter_alg(U_guess, *args)
    U_guess = U_new
    N_array.append(N)
    N0_array.append(N0)
    U_array.append(U_guess / eV)

print(f'N0 = {N0:.3e}, N = {N:.3e}')
print(f"U = {U_new / eV:.6f} eV")

fig, ax = plt.subplots(1, 3, figsize=(15, 5))
ax[0].plot(N_array[1:], '.')
ax[0].set_title("N")    
ax[1].plot(N0_array[1:], '.')
ax[1].set_title("N0")
ax[2].plot(U_array[1:], '.')
ax[2].set_title("U")
plt.show()
