from vector_utils import *

def load_U_iter(filename: str, T, V_DS, V_G, U_guess=-0.3*eV, n_iter=2000):
    if not os.path.exists(filename+'.npy'):
        U_iter = iter_U(U_guess, T, V_DS, V_G, n_iter=n_iter, delta=0.01)
        np.save(filename, U_iter)
    U_iter = np.load(filename+'.npy')
    U = U_iter[int(4*n_iter/5):].mean(axis=0)[None, :]   # (1, 2, 100, 5)
    return U_iter, U


# Parte a)
n_iter = 300                 # Es suficiente si delta = 0.01
U_guess = -0.3*eV

U_iter, U = load_U_iter('U_iter_a', T_arr, V_DS_arr, V_G_arr, U_guess=U_guess, n_iter=n_iter)
I = get_I(U, T_arr, V_DS_arr)

plt.figure(figsize=(5, 5))
plt.title('$V_G$ = $0.5$ V, $V_{DS}$ = $0.5$ V')
plt.plot(U_iter[:, 0, -1, -1]/eV, label='T = 1 K')
plt.plot(U_iter[:, 1, -1, -1]/eV, label='T = 298 K')
plt.legend()
plt.grid()
plt.xlabel('Iteración')
plt.ylabel('U [eV]')
plt.tight_layout()
plt.show()
plt.savefig('fig_a_iter.png', dpi=300)

fig, ax = plt.subplots(1, 2, figsize=(10, 5), tight_layout=True, sharey=True)
for j, T in enumerate(T_arr):
    for k, VG in enumerate(V_G_arr):
        ax[j].plot(V_DS_arr, I[0, j, :, k]/1e-6, label=f'$V_G$ = {VG:.2f} V')
    ax[j].legend()
    ax[j].set_xlabel('Drain Source bias $V_{DS}$ [V]')
    ax[j].set_title(f'$T$ = {T} K')
    ax[j].grid()
ax[0].set_ylabel('Drain Source current $I_{DS}$ [$\mu$A]')

plt.show()
# fig.savefig('img/fig_a.png', dpi=300)



# Parte b)
V_DS_arr_b = np.linspace(0., .05, 100)
U_iter_b, U_b = load_U_iter('U_iter_b', T_arr, V_DS_arr_b, V_G_arr, U_guess=U_guess, n_iter=n_iter)
I_b = get_I(U_b, T_arr, V_DS_arr_b)
I_0K_b = get_I_0K(V_DS_arr_b, V_G_arr)
I_0K = get_I_0K(V_DS_arr, V_G_arr)

for i, I_, I_0, V_DS in zip([0, 1], [I_b, I], [I_0K_b, I_0K], [V_DS_arr_b, V_DS_arr]):
    fig, ax = plt.subplots(figsize=(5, 5), tight_layout=True)
    for k, VG in enumerate(V_G_arr):
        ax.plot(V_DS, I_[0, 0, :, k]/1e-6, label=f'$V_G$ = {VG:.2f} V')
    fig.gca().set_prop_cycle(None)
    for k, VG in enumerate(V_G_arr):
        ax.plot(V_DS, I_0[0, 0, :, k]/1e-6, "--")
    ax.legend(loc='upper right')
    ax.set_xlabel('Drain Source bias $V_{DS}$ [V]')
    ax.set_title('$T$ = 0 K')
    ax.grid()
    ax.set_ylabel('Drain Source current $I_{DS}$ [$\mu$A]')

    plt.show()
#     fig.savefig(f'img/fig_b_{i}.png', dpi=300)



# Parte c)
n_iter_c = 300
V_G_arr_c = np.linspace(0, 1., 100)
V_DS_arr_c = np.array([0.5])
T_arr_c = np.array([298.])

U_iter_c, U_c = load_U_iter('U_iter_c', T_arr_c, V_DS_arr_c, V_G_arr_c, U_guess=U_guess, n_iter=n_iter_c)
I_c = get_I(U_c, T_arr_c, V_DS_arr_c)
g_m = np.diff(np.log10(I_c[0, 0, 0, :])) / np.diff(V_G_arr_c)

fig, ax = plt.subplots(figsize=(5, 5), tight_layout=True)
ax.plot(V_G_arr_c, I_c[0, 0, 0, :], label=f'$V_{{DS}}$ = {V_DS_arr_c[0]:.2f} V')
# plot the derivative
ax.plot(V_G_arr_c[1:], np.log10(g_m), label=f'$g_m$')
ax.legend()
ax.set_xlabel('Gate bias $V_G$ [V]')
ax.set_ylabel('Drain Source current $I_{DS}$ [A]')
ax.set_title(f'$T$ = {T_arr_c[0]} K')
ax.grid()
#ax.set_yscale('log')
# ax.text(0.05, 0.01, f'$g_{{m}}$ = {:.2}')
# fig.savefig('img/fig_c.png', dpi=300)

plt.show()

breakpoint()
