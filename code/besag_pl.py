# Besag pseudolikelihood estimate of the coupling in 2D SU(2) lattice gauge theory.
# Link conditional: a0 has density sqrt(1-a0^2) exp(beta k a0), mean I_2(beta k)/I_1(beta k).
import sys, numpy as np
from scipy.special import ive
from scipy.optimize import brentq
import vl_core as V
V.L = 12; V.rng = np.random.default_rng(7)
def stats(U):
    ks, a0s = [], []
    for x in range(V.L):
        for y in range(V.L):
            for mu in range(2):
                A = V.staple(U, x, y, mu); k = np.linalg.norm(A)
                ks.append(k); a0s.append(V.qmul(U[x, y, mu], A/k)[0])
    return np.array(ks), np.array(a0s)
def mple(ks, a0s):
    score = lambda b: np.sum(ks*(a0s - ive(2, b*ks)/ive(1, b*ks)))
    return brentq(score, 0.05, 50)
out = []
for beta in [1.0, 2.0, 4.0, 8.0]:
    U = V.haar((V.L, V.L, 2)); U = V.slice_heatbath(U, beta, 60)
    est = []
    for rep in range(8):
        U = V.slice_heatbath(U, beta, 10); est.append(mple(*stats(U)))
    out.append((beta, np.mean(est), np.std(est)/np.sqrt(len(est))))
    sd1 = np.std(est, ddof=1)   # spread of a single configuration's estimate
    print("beta %.1f  MPLE %.3f +- %.3f   single configuration sd %.3f (%.1f%% of beta), bias %+.3f (%+.1f se)" % (out[-1]+(sd1, 100*sd1/beta, out[-1][1]-beta, (out[-1][1]-beta)/out[-1][2])), flush=True)
import os; os.makedirs("../results", exist_ok=True); np.save("../results/besag_pl.npy", np.array(out))
