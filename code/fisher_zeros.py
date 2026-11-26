# Fisher zeros in the complex coupling of two dimensional SU(2) lattice gauge theory: a discrete system in Kac's sense.
# One plaquette (or a disc of A plaquettes): Z = (2 I_1(beta)/beta)^A, zeros at beta = +- i j_{1,n}, on the imaginary axis.
# Torus of A plaquettes: Z = sum_n (2 I_n(beta)/beta)^A, the sum over the characters of SU(2).
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.special import iv, jn_zeros
from mpmath import mp, besseli, findroot, mpc
mp.dps = 30
def F(b, A, NT=60): return sum((2*besseli(n, b)/b)**A for n in range(1, NT))
def zeros(A, re=np.linspace(-8, 8, 161), im=np.linspace(0.2, 16, 159)):
    X, Y = np.meshgrid(re, im); B = X+1j*Y; n = np.arange(1, 60)[:, None, None]
    V = np.abs(((2*iv(n, B[None])/B[None])**A).sum(0))/np.exp(A*np.abs(X))
    found = []
    for i in range(1, V.shape[0]-1):
        for j in range(1, V.shape[1]-1):
            if V[i, j] == V[i-1:i+2, j-1:j+2].min():
                try:
                    z = findroot(lambda b: F(b, A), mpc(X[i, j], Y[i, j]), tol=1e-18, maxsteps=60)
                    z = complex(z)
                    if 0.05 < z.imag < 16.5 and abs(z.real) < 8.5 and all(abs(z-w) > 1e-6 for w in found): found.append(z)
                except Exception: pass
    return np.array(found)
fig, ax = plt.subplots(figsize=(5.6, 3.6)); j = jn_zeros(1, 5)
ax.plot(np.zeros(5), j, "k*", ms=9, label="disc: $i\\,j_{1,n}$")
for A, c in [(2, "#1f4e79"), (4, "#c0392b"), (9, "#16a085")]:
    z = zeros(A); print("torus A=%d: %d zeros in the window; max |Re| %.3f; nearest to the real axis: %s" % (A, len(z), np.abs(z.real).max() if len(z) else 0, np.round(z[np.argmin(z.imag)], 3) if len(z) else None))
    ax.plot(z.real, z.imag, "o", ms=3, color=c, label=f"torus, $A={A}$")
ax.axvline(0, color="k", lw=0.5); ax.set_xlabel(r"Re $\beta$"); ax.set_ylabel(r"Im $\beta$"); ax.legend(frameon=False, fontsize=7)
fig.tight_layout(); fig.savefig("figures/fisher_zeros.pdf")
