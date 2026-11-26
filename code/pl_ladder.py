# Particle learning of the decrement of the ladder from samples of plaquettes at each rung.
# State x_m = coupling at rung m; x_{m+1} = x_m - delta + sqrt(q) eps; y_m = mean of n plaquettes at rung m.
# Parameters (delta, q) have a normal inverse gamma posterior with sufficient statistics (count, sum d, sum d^2), d = x_m - x_{m+1}.
import sys, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 0   # default seed 0 reproduces the original run
rng = np.random.default_rng(8+SEED)
M, NMAX = 3000, 300; xg, wg = np.polynomial.legendre.leggauss(M); th = (xg+1)*np.pi/2; wg = wg*np.pi/2
haar = (2/np.pi)*np.sin(th)**2*wg; nn = np.arange(1, NMAX+1); CHI = np.sin(np.outer(nn, th))/np.sin(th)
m_of = lambda x: np.exp(-1.5/x); v_of = lambda x: (1+3*np.exp(-4/x))/4-np.exp(-3/x); dm = lambda x: np.exp(-1.5/x)*1.5/x**2
# data: exact unprojected flow from the heat kernel law at beta = 8, n plaquettes sampled at each rung
u = np.exp(-(2/8.0)*(nn**2-1)/4); nobs, T = 400, 22; ys, truth = [], []
for m in range(T):
    f = np.clip((nn*u)@CHI, 1e-300, None); p = f*haar; p /= p.sum()
    ys.append(np.cos(rng.choice(th, nobs, p=p)).mean()); truth.append(-1.5/np.log(u[1]))
    g = f**4; c = (CHI@(g*haar))/nn; u = (c/c[0])**4
ys, truth = np.array(ys), np.array(truth); print("true decrements: mean %.4f, first %.4f, last %.4f" % (-np.diff(truth).mean(), truth[0]-truth[1], truth[-2]-truth[-1]))
N = 20000; x = rng.uniform(4, 12, N); w = np.exp(-0.5*(ys[0]-m_of(x))**2/(v_of(x)/nobs)); x = x[rng.choice(N, N, p=w/w.sum())]
# prior: delta | q ~ N(0.5, q/k0), q ~ IG(a0, b0)
k0, mu0, a0, b0 = 0.05, 0.5, 2.0, 0.02; cnt = np.zeros(N); s1 = np.zeros(N); s2 = np.zeros(N); post = []
def draw(cnt, s1, s2):
    k = k0+cnt; mu = (k0*mu0+s1)/k; a = a0+cnt/2; b = b0+0.5*(s2+k0*mu0**2-k*mu**2)
    q = b/rng.gamma(a); return mu+np.sqrt(q/k)*rng.standard_normal(mu.size), q
delta, q = draw(cnt, s1, s2)
for t in range(1, T):
    mu = x-delta; R = v_of(mu)/nobs; H = dm(mu); S = R+H*H*q
    lw = -0.5*(ys[t]-m_of(mu))**2/S-0.5*np.log(S); w = np.exp(lw-lw.max()); idx = rng.choice(N, N, p=w/w.sum())     # resample by the predictive
    x, mu, R, H, S, q_, cnt, s1, s2 = x[idx], mu[idx], R[idx], H[idx], S[idx], q[idx], cnt[idx], s1[idx], s2[idx]
    K = q_*H/S; xn = mu+K*(ys[t]-m_of(mu))+np.sqrt(np.clip(q_-K*H*q_, 1e-12, None))*rng.standard_normal(N)        # propagate
    d = x-xn; cnt, s1, s2 = cnt+1, s1+d, s2+d*d; x = xn; delta, q = draw(cnt, s1, s2)                              # update, redraw
    post.append((np.quantile(delta, [0.05, 0.5, 0.95]), np.quantile(x, [0.05, 0.5, 0.95]), len(np.unique(idx))/N))
D = np.array([p[0] for p in post]); X = np.array([p[1] for p in post])
print("posterior of delta after %d rungs: median %.4f, 90%% interval [%.4f, %.4f]; unique ancestors %.2f" % (T-1, D[-1, 1], D[-1, 0], D[-1, 2], post[-1][2]))
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0)); r = np.arange(1, T)
ax[0].fill_between(r, X[:, 0], X[:, 2], color="#aed6f1"); ax[0].plot(r, X[:, 1], color="#1f4e79", label="filtered coupling"); ax[0].plot(np.arange(T), truth, "k--", lw=1, label="exact flow")
ax[0].set_xlabel("rung"); ax[0].set_ylabel(r"$\beta_m$"); ax[0].legend(frameon=False, fontsize=7)
ax[1].fill_between(r, D[:, 0], D[:, 2], color="#f5b7b1"); ax[1].plot(r, D[:, 1], color="#c0392b", label=r"posterior of $\delta$"); ax[1].axhline(-np.diff(truth).mean(), color="k", ls="--", lw=1, label="mean true decrement")
ax[1].set_ylim(0, 0.8); ax[1].set_xlabel("rung"); ax[1].set_ylabel("decrement"); ax[1].legend(frameon=False, fontsize=7)
fig.tight_layout()
if SEED == 0: fig.savefig("../figures/pl_ladder.pdf")
