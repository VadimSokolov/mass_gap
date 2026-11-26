# The mean count as an observable: E[1/(1 - cos theta)] - 2, and its expansion 2 S_1 - 3/2 + 3 S_2/S_1 on the class L,
# with S_1 = sum 1/(gamma_m + 1), S_2 = sum 1/(gamma_m + 1)^2.  Then the first rung against -13/16 + (15/16) r, r = 1 - 2 S_2/S_1.
import numpy as np
from scipy.special import gammaln
from count_lyapunov import MM, U, logZ, rung, j
th = (np.arange(400000)+0.5)*np.pi/400000; x = np.cos(th); haar = np.sin(th)**2
def observable(logrho): w = np.exp(logrho-logrho.max())*haar; return (w/(1-x)).sum()/w.sum()-2
def count_law(gam):
    nu = np.array([1.0])
    for g in gam: nu = np.convolve(nu, [1-1/g, 1/g])
    l = np.full(MM+1, -np.inf); ok = nu > 0; l[:nu.size][ok] = np.log(nu[ok])+logZ[:nu.size][ok]; p = np.exp(l-l.max()); return p/p.sum()
if __name__ == "__main__":
    rng = np.random.default_rng(2)
    print("f_j:", [round(float(observable(jj*np.log1p(x))), 6) for jj in (1, 5, 40)], " Wilson kappa = 4, 16:", [round(float(observable(k*x)), 4) for k in (4, 16)], "against 2 kappa - 3/2")
    print("%-34s %9s %9s %9s %8s %9s %9s" % ("law in L", "count", "observ.", "expans.", "r", "rung", "law"))
    for name, gam in [("40 zeros at gamma = 1", np.ones(40)), ("60 zeros at gamma = 2", np.full(60, 2.0)), ("120 zeros at gamma = 5", np.full(120, 5.0)), ("300 zeros at gamma = 20", np.full(300, 20.0)),
                      ("80 zeros uniform on [1, 6]", rng.uniform(1, 6, 80)), ("150 zeros, gamma = 1 + exponential", 1+rng.exponential(3, 150))]:
        pi = count_law(gam); S1 = (1/(gam+1)).sum(); S2 = (1/(gam+1)**2).sum(); r = 1-2*S2/S1; m0, mb, m1, u2 = rung(pi)
        print("%-34s %9.4f %9.4f %9.4f %8.4f %+9.4f %+9.4f" % (name, m0, observable(np.log1p(x[:, None]/gam[None, :]).sum(1)), 2*S1-1.5+3*S2/S1, r, m1-m0, -13/16+15/16*r))
