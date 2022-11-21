from vector_utils import *



# Parte a)
n_iter = 2000                # Es suficiente si delta = 0.01
U_guess = -0.3*eV

if not os.path.exists('U_iter.npy'):
    U_iter = iter_U(U_guess, T_array, V_DS_array, V_G_array, n_iter=n_iter, delta=0.01)
    np.save('U_iter', U_iter)
U_iter = np.load('U_iter.npy')
U = U_iter[int(4*n_iter/5):].mean(axis=0)[None, :]   # (1, 2, 100, 5)
I = get_I(U, T_array, V_DS_array)

# plot_U(U, U_iter)

fig, ax = plt.subplots(1, 2, figsize=(10, 5), tight_layout=True, sharey=True)
for j, T in enumerate(T_array):
    for k, VG in enumerate(V_G_array):
        ax[j].plot(V_DS_array, I[0, j, :, k]/1e-6, label=f'$V_G$ = {VG:.2f} V')
    ax[j].legend()
    ax[j].set_xlabel('Drain Source bias $V_{DS}$ [V]')
    ax[j].set_title(f'$T$ = {T} K')
    ax[j].grid()
ax[0].set_ylabel('Drain Source current $I_{DS}$ [$\mu$A]')

plt.show()
fig.savefig('img/fig_a.pdf')



# Parte b)
I_0K = get_I_0K(V_DS_array, V_G_array)

fig, ax = plt.subplots(figsize=(5, 5), tight_layout=True)
for k, VG in enumerate(V_G_array):
    ax.plot(V_DS_array, I[0, 0, :, k]/1e-6, label=f'$V_G$ = {VG:.2f} V')
fig.gca().set_prop_cycle(None)
for k, VG in enumerate(V_G_array):
    ax.plot(V_DS_array, I_0K[0, 0, :, k]/1e-6, "--")
ax.legend(loc='upper right')
ax.set_xlabel('Drain Source bias $V_{DS}$ [V]')
ax.set_title('$T$ = 0 K')
ax.grid()
ax.set_ylabel('Drain Source current $I_{DS}$ [$\mu$A]')

plt.show()
fig.savefig('img/fig_b.pdf')



# Parte c)
n_iter2 = 200
V_G_arr2 = np.linspace(0, 1., 100)
V_DS_arr2 = np.array([0.5])
T_arr2 = np.array([298])

if not os.path.exists('U2_iter.npy'):
    U2_iter = iter_U(U_guess, T_arr2, V_DS_arr2, V_G_arr2, n_iter=n_iter2, delta=0.01)
    np.save('U2_iter', U2_iter)
U2_iter = np.load('U2_iter.npy')
U2 = U2_iter[int(4*n_iter2/5):].mean(axis=0)[None, :]   # (1, 2, 100, 5)
I2 = get_I(U2, T_arr2, V_DS_arr2)
g_m = np.diff(I2[0, 0, 0, :]) / np.diff(V_G_arr2)

fig, ax = plt.subplots(figsize=(5, 5), tight_layout=True)
ax.plot(V_G_arr2, I2[0, 0, 0, :], label=f'$V_{{DS}}$ = {V_DS_arr2[0]:.2f} V')
# plot the derivative
# ax.plot(V_G_arr2[1:], g_m/1e-6, label=f'$g_m$')
ax.legend()
ax.set_xlabel('Gate bias $V_G$ [V]')
ax.set_ylabel('Drain Source current $I_{DS}$ [A]')
ax.set_title(f'$T$ = {T_arr2[0]} K')
ax.grid()
ax.set_yscale('log')
# ax.text(0.05, 0.01, f'$g_{{m}}$ = {:.2}')
fig.savefig('img/fig_c.pdf')

plt.show()

breakpoint()
