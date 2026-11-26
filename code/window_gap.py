# One sided windows and the gap, on a lattice of spacing a in Euclidean time.
# For C(tau) = sum_n w_n exp(-E_n tau) with w_n >= 0, the negative binomial moments
#   A_k(q) = sum_{j>=0} C(ja) binom(j+k, k) q^j = sum_n w_n / (1 - q exp(-E_n a))^(k+1)
# give M_k = -log((1 - A_{k-1}/A_k)/q)/a, which decreases to the gap as k grows: a one sided, monotone estimator (the lattice form of the
# Post and Widder kernel of the window theorem). Test: one plaquette at lambda = 1/2, operator cos(theta), exact correlator plus noise,
# against the usual effective mass.
import numpy as np
from scipy.special import gammaln
from quantum_plaquette import ham
lam, N = 0.5, 80; e, v = np.linalg.eigh(ham(lam, N)); X = np.diag(np.full(N-1, 0.5), 1)+np.diag(np.full(N-1, 0.5), -1)
w = (v.T@X@v)[0, 1:]**2; E = e[1:]-e[0]; gap = E[0]; a = 0.25; J = np.arange(0, 161); tau = a*J; C = (w[None, :]*np.exp(-np.outer(tau, E))).sum(1)
def nb(Cn, k, sigma):
    q = np.exp(-sigma*a); A = lambda m: (Cn*np.exp(gammaln(J+m+1)-gammaln(J+1)-gammaln(m+1)+J*np.log(q))).sum(); r = 1-A(k-1)/A(k); return -np.log(r/q)/a if r > 0 else np.nan
def eff(Cn, t, lag=1.0): i, j = int(round(t/a)), int(round((t+lag)/a)); return -np.log(Cn[j]/Cn[i])/lag if Cn[i] > 0 and Cn[j] > 0 else np.nan
if __name__ == "__main__":
    print("gap %.5f; weights of the first three states %s" % (gap, np.round(w[:3]/w.sum(), 5)))
    print("exact, sigma = 3: k: M_k   ", ["%d: %.5f" % (k, nb(C, k, 3.0)) for k in (1, 2, 4, 8, 16, 32)])
    print("exact,            t: m_eff ", ["%g: %.5f" % (t, eff(C, t)) for t in (0.5, 1, 2, 3, 4)])
    rng = np.random.default_rng(6); noise = 1e-4*C[0]; R = 2000; draw = lambda: C+noise*rng.standard_normal(C.size)
    for sigma, k in ((3.0, 6), (3.0, 12), (3.0, 20), (3.0, 32)):
        s = np.array([nb(draw(), k, sigma) for _ in range(R)]); print("negative binomial, sigma = %g, k = %2d: exact %.5f   with noise %.5f +- %.5f" % (sigma, k, nb(C, k, sigma), np.nanmean(s), np.nanstd(s)))
    for t in (1.0, 2.0, 3.0, 4.0):
        s = np.array([eff(draw(), t) for _ in range(R)]); print("effective mass, t = %g: exact %.5f   with noise %.5f +- %.5f" % (t, eff(C, t), np.nanmean(s), np.nanstd(s)))
