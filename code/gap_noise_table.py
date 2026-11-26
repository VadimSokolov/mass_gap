"""Support constrained profile estimate of the gap at several noise levels (toy spectrum)."""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import nnls
atoms = [(1.0, 1.0), (1.5, 0.6), (1.9, 0.3)]
rho = lambda E: 0.4*np.sqrt(max(E-2, 0))*np.exp(-0.3*(E-2))
tau = np.arange(0.1, 6.01, 0.1); Eg = np.linspace(0.2, 10, 246)
C = np.array([sum(c*np.exp(-t*E) for E, c in atoms)+quad(lambda E: rho(E)*np.exp(-t*E), 2, 60)[0] for t in tau])
K = np.exp(-np.outer(tau, Eg)); cands = np.linspace(0.5, 1.4, 91)
for rel in [0.001, 0.002, 0.005, 0.01, 0.02]:
    est = []
    for rep in range(20):
        rng = np.random.default_rng(100+rep); sd = rel*C; y = C + sd*rng.standard_normal(len(C))
        chi = np.array([nnls((K/sd[:, None])[:, Eg >= d], y/sd)[1]**2 for d in cands])
        est.append(cands[np.argmax(chi-chi.min() > 4.0)])
    est = np.array(est); print("noise %.1f%%: mean %.3f  sd %.3f  min %.3f max %.3f" % (100*rel, est.mean(), est.std(), est.min(), est.max()))
