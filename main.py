from utils import *

U_guess = -.311*eV                                  # correcto -0.311 *eV
V_DS = V_DS_array.mean()
U_ana = deme_un_U_mi_Rey(V_G_array)                 # Válido para T = 0 K

I_data = np.empty((len(T_array), len(V_G_array), len(V_DS_array)))

fig, ax = plt.subplots(1, 2, figsize=(10, 5), tight_layout=True, sharey=True)
for j, T in enumerate(T_array):
    for k, VG in enumerate(V_G_array):
        U = get_U_iterative(U_guess, T, V_DS, VG)[-500:].mean()
        I_array = np.array([calculate_I(U, T, VDS) for VDS in V_DS_array])
        I_data[j, k, :] = I_array
        ax[j].plot(V_DS_array, I_array/1e-6, label=f'$V_G$ = {VG:.2f} V')
        I_ana = np.array([calculate_I(U_ana[k], T, VDS) for VDS in V_DS_array])
        ax[j].plot(V_DS_array, I_ana/1e-6, '--', c='k', label=f'$V_G$ = {VG:.2f} V (analytical)')
    ax[j].legend()
    ax[j].set_xlabel('$V_{DS}$ [V]')
    ax[j].set_title(f'T = {T} K')
    ax[j].grid()
ax[0].set_ylabel('$I_{DS}$ [$\mu$A]')

plt.show()
fig.savefig('img/fig_1.pdf')


