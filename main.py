from vector_utils import *

n_iter = 2000                # Es suficiente si delta = 0.01
U_guess = -0.3*eV

if not os.path.exists('U_iter.npy'):
    U = iter_U(U_guess, T_array, V_DS_array, V_G_array, n_iter=n_iter, delta=0.01)
    np.save('U_iter', U_iter)
U_iter = np.load('U_iter.npy')
U = U_iter[int(4*n_iter/5):].mean(axis=0)[None, :]   # (1, 2, 100, 5)
I = get_I(U, T_array, V_DS_array)

plot_U(U, U_iter)

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
fig.savefig('img/fig1.pdf')

breakpoint()

