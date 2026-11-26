# Is the mean count a Lyapunov function for the exact ladder (k = 4 copies, fourfold convolution) on the whole family M?
# Laws are handled through characters: rho = sum_j pi(j) f_j has u_n = sum_j pi(j) u_n(j+1) and mean count 2 sum_n n u_n - 2.
# Bond moving acts on nu_0 (pi with the Catalan weights removed) by fourfold convolution; decimation raises u_n to the fourth power.
import numpy as np
from scipy.special import gammaln
MM = 400; J = np.arange(1, MM+2); n = np.arange(1, MM+2)
with np.errstate(invalid="ignore"): LU = gammaln(J)[None, :]+gammaln(J+2)[None, :]-gammaln(np.maximum(J[None, :]-n[:, None], 0)+1)-gammaln(J[None, :]+n[:, None]+1)
U = np.where(n[:, None] <= J[None, :], np.exp(LU), 0.0); j = np.arange(MM+1); logZ = gammaln(2*j+3)-2*gammaln(j+2)-np.log(j+2)-j*np.log(2)
mean_of_u = lambda u: 2*(n*u).sum()-2
def rung(pi, k=4, power=4):
    lnu = np.where(pi > 0, np.log(np.where(pi > 0, pi, 1))-logZ, -np.inf); lnu -= lnu.max(); nu = np.exp(lnu); nu /= nu.sum(); out = nu.copy()
    for _ in range(k-1): out = np.convolve(out, nu)[:MM+1]
    lp = np.where(out > 0, np.log(np.where(out > 0, out, 1))+logZ, -np.inf); p = np.exp(lp-lp.max()); p /= p.sum(); u = U@p; return mean_of_u(U@pi), (p*j).sum(), mean_of_u(u**power), (u**power)[1]
def pad(v): out = np.zeros(MM+1); out[:len(v)] = v; return out/np.sum(v)
if __name__ == "__main__":
    rng = np.random.default_rng(1); worst = (-9, None); rows = []
    fam = {"point mass at %d" % a: pad([0]*a+[1]) for a in (1, 2, 5, 20, 60)}
    fam.update({"Poisson nu_0, kappa = %g (Wilson)" % k: None for k in (0.5, 1, 2, 4, 8, 16)})
    for name in list(fam):
        if fam[name] is None:
            k = float(name.split("= ")[1].split(" ")[0]); l = j*np.log(k)-gammaln(j+1)+logZ; fam[name] = pad(np.exp(l-l.max()))
    fam.update({"two points 0 and %d, weight %.2f on 0" % (a, w): pad([w]+[0]*(a-1)+[1-w]) for a in (3, 20, 60) for w in (0.5, 0.9, 0.99)})
    fam.update({"geometric, ratio %.2f" % r: pad(r**np.arange(90)) for r in (0.5, 0.9, 0.97)})
    for name, pi in fam.items():
        m0, mb, m1, u2 = rung(pi); rows.append((name, m0, mb, m1)); print("%-36s mean count %8.3f  after bond moving %8.3f  after the rung %8.3f  change %+8.4f  u_2 change %+.4f" % (name, m0, mb, m1, m1-m0, u2-(pi*j/(j+3)).sum()))
    for trial in range(4000):
        size = rng.integers(2, 60); pi = pad(rng.dirichlet(np.full(size, rng.choice([0.05, 0.3, 1.0]))))
        m0, mb, m1, u2 = rung(pi)
        if m1-m0 > worst[0]: worst = (m1-m0, pi[:size].round(3))
    print("largest change of the mean count over 4000 random laws: %+.5f" % worst[0])
    u2c = []
    for trial in range(4000):
        size = rng.integers(2, 60); pi = pad(rng.dirichlet(np.full(size, rng.choice([0.05, 0.3, 1.0])))); m0, mb, m1, u2 = rung(pi); u2c.append(u2-(pi*j/(j+3)).sum())
    print("largest change of u_2 over 4000 random laws: %+.6f" % max(u2c))
    # log concave nu_0: sums of independent Bernoulli variables, and exponentials of random concave sequences
    res = {"bernoulli sums": [], "concave": []}
    for trial in range(3000):
        if trial % 2 == 0:
            nu = np.array([1.0])
            for q in rng.uniform(0, 1, rng.integers(1, 40))**rng.choice([1, 3]): nu = np.convolve(nu, [1-q, q])
            kind = "bernoulli sums"
        else:
            L = rng.integers(2, 50); nu = np.exp(np.cumsum(np.concatenate([[0], np.sort(rng.normal(0, rng.choice([0.3, 1, 3]), L))[::-1]]))); kind = "concave"
        pi = pad(nu*np.exp(logZ[:nu.size]-logZ[:nu.size].max())); m0, mb, m1, u2 = rung(pi); res[kind].append((m1-m0, u2-(pi*j/(j+3)).sum(), m0))
    for kind, v in res.items(): v = np.array(v); print("log concave nu_0 (%s): largest change of mean count %+.4f (at mean %.1f), largest change of u_2 %+.2e, over %d laws" % (kind, v[:, 0].max(), v[np.argmax(v[:, 0]), 2], v[:, 1].max(), len(v)))
