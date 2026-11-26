# Run the ladder and perturbation theory in Parisi's effective coupling beta_E = 3/(4(1 - <P>)),
# with <P> from our own 8^4 heat bath runs (su2_4d.py), against SU(2) string tension data.
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from mk_verify import make_R
beta = np.array([2.2, 2.3, 2.4, 2.5, 2.635]); asig = np.array([0.4690, 0.3690, 0.2660, 0.1871, 0.1207])
P = np.array([0.569972, 0.603046, 0.629759, 0.652635, 0.676289]); dP = np.array([3.75, 3.51, 2.93, 2.00, 1.23])*1e-4
bE = 3/(4*(1-P)); print("beta_E", np.round(bE, 4))
R = make_R(3000, 600)
N = 2; b0 = 11*N/(48*np.pi**2); b1 = 34*N**2/(3*(16*np.pi**2)**2)
def logA2(b): g2 = 2*N/b; return -b1/(2*b0**2)*np.log(b0*g2) - 1/(2*b0*g2)
# ladder: log a decreases by log 2 per decrement delta(beta) = beta - R(beta); integrate d log a = -log2 dbeta/delta
grid = np.linspace(1.0, 3.5, 251); rate = np.array([np.log(2)/(g-R(g)) for g in grid])
cum = np.concatenate([[0], np.cumsum(0.5*(rate[1:]+rate[:-1])*np.diff(grid))])
logA_mk = lambda b: -np.interp(b, grid, cum)
for name, x in [("bare", beta), ("effective", bE)]:
    s_data = np.polyfit(x, np.log(asig), 1)[0]
    s_mk = np.polyfit(x, logA_mk(x), 1)[0]; s_2l = np.polyfit(x, logA2(x), 1)[0]
    r_mk = np.std(np.log(asig) - logA_mk(x) - np.mean(np.log(asig) - logA_mk(x)))
    r_2l = np.std(np.log(asig) - logA2(x) - np.mean(np.log(asig) - logA2(x)))
    print(f"{name:9s} data {s_data:.3f}  ladder {s_mk:.3f}  two loop {s_2l:.3f}  resid ladder {r_mk:.4f} two loop {r_2l:.4f}")
print("decrements", [round(g-R(g), 4) for g in [1.5, 1.7, 1.9, 2.1]])
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.3), sharey=True)
for a, x, t in [(ax[0], beta, r'bare $\beta$'), (ax[1], bE, r'effective $\beta_E$')]:
    xx = np.linspace(x.min()-0.03, x.max()+0.03, 60)
    a.semilogy(x, asig, 'o', color='#c0392b', label='lattice data')
    for f, c, l in [(logA_mk, '#1f4e79', 'ladder'), (logA2, '#e67e22', 'two loop')]:
        a.semilogy(xx, np.exp(f(xx) - f(x[3]))*asig[3], color=c, label=l)
    a.set_xlabel(t)
ax[0].set_ylabel(r'$a\sqrt{\sigma}$'); ax[0].legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig('../figures/effective_coupling.pdf')
