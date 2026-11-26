# Does the SU(3) heat kernel have monotone likelihood ratio in time along rays from the identity in the Weyl alcove?
# K_t(theta) = sum_{p,q} dim(p,q) exp(-t C2(p,q)) Re chi_{pq}(theta), theta = r h on the maximal torus.
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
P = 60; reps = [(p, q) for p in range(P+1) for q in range(P+1-p)]
C2 = np.array([(p*p+q*q+p*q+3*p+3*q)/3 for p, q in reps]); dim = np.array([(p+1)*(q+1)*(p+q+2)/2 for p, q in reps])
def chars(th):                      # th: (3, n) angles summing to zero; returns (nreps, n) real parts of characters
    E = lambda l, k: np.exp(1j*l*th[k])
    det = lambda l: (E(l[0],0)*(E(l[1],1)*E(l[2],2)-E(l[1],2)*E(l[2],1)) - E(l[0],1)*(E(l[1],0)*E(l[2],2)-E(l[1],2)*E(l[2],0)) + E(l[0],2)*(E(l[1],0)*E(l[2],1)-E(l[1],1)*E(l[2],0)))
    dn = det([2, 1, 0]); return np.array([np.real(det([p+q+2, q+1, 0])/dn) for p, q in reps])
def logK(t, CH): return np.log(np.clip((dim*np.exp(-t*C2))@CH, 1e-300, None))
dirs = {"(1, 0, -1)": (1, 0, -1), "(1, 0.5, -1.5)": (1, .5, -1.5), "(1, 0.9, -1.9)": (1, .9, -1.9), "(1, -0.3, -0.7)": (1, -.3, -.7)}
ts = np.geomspace(0.12, 6, 70); fig, ax = plt.subplots(1, 4, figsize=(7.6, 2.6), sharey=True); worst = {}
for a, (name, h) in zip(ax, dirs.items()):
    h = np.array(h, float); rmax = 2*np.pi/(h.max()-h.min()); r = np.linspace(0.02, 0.98, 120)*rmax
    CH = chars(np.outer(h, r)); CHp = chars(np.outer(h, r+1e-4)); CHm = chars(np.outer(h, r-1e-4))
    Z = np.array([((logK(t*(1+1e-4), CHp)-logK(t*(1-1e-4), CHp))-(logK(t*(1+1e-4), CHm)-logK(t*(1-1e-4), CHm)))/(2e-4*t*2e-4) for t in ts])
    ok = np.array([(dim*np.exp(-t*C2))@CH for t in ts]) > 1e-6; Z = np.where(ok, Z, np.nan); worst[name] = (np.nanmin(Z), ok.mean()); pc = a.pcolormesh(r/rmax, ts, np.sign(Z)*np.log10(1+np.abs(Z)), cmap="viridis", vmin=-1, vmax=3, shading="auto")
    a.set_yscale("log"); a.set_title("ray "+name, fontsize=7); a.set_xlabel("$r/r_{\\max}$")
ax[0].set_ylabel("$t$"); fig.colorbar(pc, ax=ax, label="signed log of the mixed partial", fraction=0.02)
fig.savefig("figures/su3_mlr.pdf", bbox_inches="tight")
for k, v in worst.items(): print("ray %-16s min mixed partial %.4g on the %.0f%% of the grid where K > 1e-6" % (k, v[0], 100*v[1]))
