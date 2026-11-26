"""Figures for Thorin measures and the Yang Mills mass gap. Polson and Sokolov, Nov 2026."""
import numpy as np, mpmath as mp, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import nnls
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
out = "../figures/"
Delta = 1.0
atoms = [(1.0, 1.0), (1.5, 0.6), (1.9, 0.3)]
def rho(E): return np.where(E > 2*Delta, 0.4*np.sqrt(np.clip(E-2*Delta, 0, None))*np.exp(-0.3*(E-2*Delta)), 0.0)
def T(alpha, z):
    y = np.sqrt(z/2); t = np.sqrt(2*z)
    s = sum(c*np.exp(-alpha*E)*2*np.sin(E*y)**2 for E, c in atoms)
    s += quad(lambda E: rho(E)*np.exp(-alpha*E)*2*np.sin(E*y)**2, 2*Delta, 60, limit=500)[0]
    return s/(np.pi*t)

# Figure 1: Thorin densities and the gap
al = np.linspace(0.5, 30, 60)
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
for z, col in [(0.8, "C0"), (3.0, "C1"), (12.0, "C2")]:
    lt = np.array([np.log(T(a, z)) for a in al])
    ax[0].plot(al, lt, col, label="z = %g" % z)
    ax[1].plot(al[1:], -np.diff(lt)/np.diff(al), col, label="z = %g" % z)
ax[0].plot(al, -Delta*al - 1.0, "k:", label="slope $-\\Delta$")
ax[0].set_xlabel("$\\alpha$ (Euclidean time)"); ax[0].set_ylabel("$\\log T_\\alpha(z)$"); ax[0].legend(frameon=False)
ax[1].axhline(Delta, color="k", ls=":"); ax[1].set_ylim(0.9, 1.6)
ax[1].set_xlabel("$\\alpha$"); ax[1].set_ylabel("$-\\partial_\\alpha \\log T_\\alpha(z)$"); ax[1].legend(frameon=False)
fig.tight_layout(); fig.savefig(out+"thorin_gap.pdf"); plt.close()

# Figure 2: Riemann kernel from the SU(2) sphere partition function
def Phi(u, N=40): return float(sum((2*mp.pi**2*n**4*mp.e**(4.5*u)-3*mp.pi*n**2*mp.e**(2.5*u))*mp.e**(-mp.pi*n**2*mp.e**(2*u)) for n in range(1, N)))
def F(A, N=40): return float(sum(n**2*mp.e**(-A*n**2/8) for n in range(1, N)))
def Fp(A, N=40): return float(sum(-n**4/8*mp.e**(-A*n**2/8) for n in range(1, N)))
uu = np.linspace(-0.6, 1.0, 200)
ph = [Phi(u) for u in uu]
ym = [2*np.pi**2*np.exp(4.5*u)*(-8*Fp(8*np.pi*np.exp(2*u))) - 3*np.pi*np.exp(2.5*u)*F(8*np.pi*np.exp(2*u)) for u in uu]
fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
ax[0].plot(uu, ph, "C0", lw=3, alpha=0.5, label="Riemann $\\Phi(u)$")
ax[0].plot(uu, ym, "k--", lw=1.2, label="from $Z_0^{SU(2)}(A)$, $A = 8\\pi e^{2u}$")
ax[0].set_xlabel("$u$"); ax[0].legend(frameon=False)
AA = np.linspace(2, 60, 200)
ax[1].semilogy(AA, [np.exp(A/8)*F(A) for A in AA], "C3")
ax[1].set_xlabel("area $A$"); ax[1].set_ylabel("$Z_0(A)$, SU(2) on the sphere")
fig.tight_layout(); fig.savefig(out+"riemann_kernel_ym.pdf"); plt.close()

# Figure 3: positive deconvolution of a noisy correlator, gap as a support constraint
rng = np.random.default_rng(7)
tau = np.arange(0.1, 6.01, 0.1)
Eg = np.linspace(0.2, 10, 246); dE = Eg[1]-Eg[0]
def Ctrue(t): return sum(c*np.exp(-t*E) for E, c in atoms) + quad(lambda E: rho(E)*np.exp(-t*E), 2*Delta, 60)[0]
y0 = np.array([Ctrue(t) for t in tau]); sd = 2e-3*y0
y = y0 + sd*rng.standard_normal(len(tau))
K = np.exp(-np.outer(tau, Eg))
Kw, yw = K/sd[:, None], y/sd
nu = np.ones(len(Eg))*1e-2
for _ in range(20000):                      # Vardi Lee EM (multiplicative, nonnegative)
    pred = K@nu
    nu *= (K.T@(np.clip(y, 1e-300, None)/pred))/K.sum(0)
cands = np.linspace(0.5, 1.4, 46); chi = []
for d in cands:
    keep = Eg >= d
    c, r = nnls(Kw[:, keep], yw); chi.append(r**2)
chi = np.array(chi)
fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
ax[0].bar(Eg, nu, width=dE, color="C0", alpha=0.7, label="EM reconstruction")
for E, c in atoms: ax[0].axvline(E, color="k", ls=":", lw=1)
ax[0].axvspan(2*Delta, 10, color="C1", alpha=0.08, label="true continuum")
ax[0].set_xlim(0, 6); ax[0].set_xlabel("$E$"); ax[0].set_ylabel("mass"); ax[0].legend(frameon=False)
ax[1].plot(cands, chi - chi.min(), "C3o-", ms=3)
ax[1].axvline(Delta, color="k", ls=":"); ax[1].set_yscale("symlog", linthresh=1); ax[1].set_ylim(0, None)
ax[1].set_xlabel("imposed gap $\\delta$ (support of $\\mu$ in $[\\delta,\\infty)$)"); ax[1].set_ylabel("$\\chi^2(\\delta) - \\min \\chi^2$")
fig.tight_layout(); fig.savefig(out+"deconvolution_gap.pdf"); plt.close()
j = np.argmax(chi - chi.min() > 4.0)
print("profile: chi2 rises by 4 at delta =", cands[j], "true gap", Delta)

# Figure 4: 2D SU(2), Thorin slope converges to the Casimir gap 3/8
def T2d(A0, z, g, N=80):
    y = np.sqrt(z/2); t = np.sqrt(2*z)
    return sum(n**(2-2*g)*np.exp(-A0*(n*n-1)/8)*2*np.sin((n*n-1)/8*y)**2 for n in range(2, N))/(np.pi*t)
A = np.linspace(1, 60, 120)
fig, ax = plt.subplots(figsize=(5.5, 3.4))
for g, col in [(0, "C0"), (1, "C2")]:
    lt = np.log([T2d(a, 5.0, g) for a in A])
    ax.plot(A[1:], -np.diff(lt)/np.diff(A), col, label="genus %d" % g)
ax.axhline(3/8, color="k", ls=":", label="gap $C_2/2 = 3/8$")
ax.set_ylim(0.3, 1.0); ax.set_xlabel("area $A$"); ax.set_ylabel("$-\\partial_A \\log T_A(z)$, $z = 5$"); ax.legend(frameon=False)
fig.tight_layout(); fig.savefig(out+"ym2d_gap.pdf"); plt.close()
print("done")
