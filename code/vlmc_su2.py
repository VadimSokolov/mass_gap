"""The vertical likelihood family (Polson and Scott) on 2D SU(2) lattice gauge theory, 3x3 torus,
target Z(beta) = E_Haar exp(-beta S), exact by characters. Methods:
 1 importance sampling from the prior          2 harmonic mean from posterior draws
 3 annealing / tempering ladder (thermodynamic integration, trapezoid)   4 slice sampling heat bath (posterior kernel)
 5 nested sampling (Skilling)                  6 energy level sampling (Wang Landau density of states)
 7 vertical likelihood reweighted sampling: draw x from g(x) proportional to p(x) w(S(x)), reweight by 1/w
Polson and Sokolov, Nov 2026."""
import numpy as np, time
from scipy.special import ive, logsumexp
rng = np.random.default_rng(5); L = 3; P = L*L
def qmul(a, b):
    a0, a1, a2, a3 = a[..., 0], a[..., 1], a[..., 2], a[..., 3]; b0, b1, b2, b3 = b[..., 0], b[..., 1], b[..., 2], b[..., 3]
    return np.stack([a0*b0-a1*b1-a2*b2-a3*b3, a0*b1+a1*b0+a2*b3-a3*b2, a0*b2-a1*b3+a2*b0+a3*b1, a0*b3+a1*b2-a2*b1+a3*b0], -1)
def conj(a): return a*np.array([1., -1, -1, -1])
def haar(shape):
    v = rng.standard_normal(shape+(4,)); return v/np.linalg.norm(v, axis=-1, keepdims=True)
def action_all(U):
    a = U[:, :, 0]; b = np.roll(U[:, :, 1], -1, axis=0); c = conj(np.roll(U[:, :, 0], -1, axis=1)); d = conj(U[:, :, 1])
    return float(np.sum(1-qmul(qmul(qmul(a, b), c), d)[..., 0]))
def staple(U, x, y, mu):
    """sum of the two staples so that the plaquettes containing link (x,y,mu) are tr(U_link * A)/2 summed"""
    if mu == 0:
        s1 = qmul(qmul(U[(x+1) % L, y, 1], conj(U[x, (y+1) % L, 0])), conj(U[x, y, 1]))
        ym = (y-1) % L
        s2 = qmul(qmul(conj(U[(x+1) % L, ym, 1]), conj(U[x, ym, 0])), U[x, ym, 1])
    else:
        xm = (x-1) % L
        s1 = qmul(qmul(conj(U[(x+1) % L, y, 0]) if False else conj(U[x, (y+1) % L, 0]), conj(U[xm, y, 1])), U[xm, y, 0])
        s2 = qmul(qmul(U[x, (y+1) % L, 0], conj(U[(x+1) % L, y, 1])), conj(U[x, y, 0]))
    return s1+s2
def logZ_exact(b):
    n = np.arange(1, 400); return -b*P + logsumexp(P*(np.log(np.maximum(2*ive(n, b)/b, 1e-300)) + b))
def meanS_exact(b, h=1e-5): return -(logZ_exact(b+h)-logZ_exact(b-h))/(2*h)
def local_S(U, x, y, mu):   # contribution of the two plaquettes containing the link
    A = staple(U, x, y, mu); return 2 - qmul(U[x, y, mu], A)[0]
def check_staple():
    U = haar((L, L, 2)); S0 = action_all(U)
    for _ in range(20):
        x, y, mu = rng.integers(L), rng.integers(L), rng.integers(2)
        old = local_S(U, x, y, mu); U2 = U.copy(); U2[x, y, mu] = haar(()); new = local_S(U2, x, y, mu)
        assert abs((action_all(U2)-S0) - (new-old)) < 1e-10, "staple check failed"
    return True
def slice_heatbath(U, beta, sweeps):
    """method 4: exact heat bath via slice sampling on a0. conditional of link: exp(beta k <W, V>) on S^3"""
    for _ in range(sweeps):
        for x in range(L):
            for y in range(L):
                for mu in range(2):
                    A = staple(U, x, y, mu); k = np.linalg.norm(A); V = conj(A)/k      # want W = U A/k
                    a0 = qmul(U[x, y, mu], A/k)[0]; kap = beta*k
                    # slice: u ~ U(0, e^{kap(a0-1)}), then a0 ~ sqrt(1-a0^2) on (a0min, 1)
                    lo = max(-1.0, a0 + np.log(rng.random())/kap)
                    while True:   # sample sqrt(1-t^2) on (lo,1) by rejection from uniform
                        t = lo + (1-lo)*rng.random()
                        if rng.random() < np.sqrt(max(1-t*t, 0)): break
                    d = rng.standard_normal(3); d /= np.linalg.norm(d)
                    W = np.concatenate([[t], np.sqrt(max(1-t*t, 0))*d])
                    U[x, y, mu] = qmul(W, conj(A/k))
    return U
def metro_constrained(U, S, Smax, nprop, eps, wfun=None):
    acc = 0
    for _ in range(nprop):
        x, y, mu = rng.integers(L), rng.integers(L), rng.integers(2)
        old = local_S(U, x, y, mu); keep = U[x, y, mu].copy()
        v = rng.standard_normal(3)*eps; th = np.linalg.norm(v)
        U[x, y, mu] = qmul(np.concatenate([[np.cos(th)], np.sin(th)*v/max(th, 1e-12)]), keep)
        Snew = S - old + local_S(U, x, y, mu)
        ok = Snew < Smax if wfun is None else np.log(rng.random()) < wfun(Snew) - wfun(S)
        if ok: S = Snew; acc += 1
        else: U[x, y, mu] = keep
    return U, S, acc/nprop
if __name__ == "__main__":
    check_staple(); out = {}; b0 = 2.0; exact = logZ_exact(b0); print("exact log Z(2) =", exact)
    # 1 prior importance sampling
    t = time.time(); Sp = np.array([action_all(haar((L, L, 2))) for _ in range(40000)])
    out['prior IS'] = logsumexp(-b0*Sp) - np.log(len(Sp)); print("1 prior IS", out['prior IS'], "%.0fs" % (time.time()-t))
    # 4 slice heat bath posterior draws at beta = 2, used by harmonic mean
    t = time.time(); U = haar((L, L, 2)); U = slice_heatbath(U, b0, 50); post = []
    for i in range(3000): U = slice_heatbath(U, b0, 1); post.append(action_all(U))
    post = np.array(post); print("4 slice heat bath: <S> =", post.mean(), " exact", meanS_exact(b0), "%.0fs" % (time.time()-t))
    # 2 harmonic mean: 1/Z = E_post[exp(beta S)]
    out['harmonic mean'] = -(logsumexp(b0*post) - np.log(len(post))); print("2 harmonic mean", out['harmonic mean'])
    # 3 tempering ladder: log Z(b0) = -int_0^b0 <S>_b db, heat bath at each rung
    t = time.time(); bs = np.linspace(0, b0, 11); ms = []
    for b in bs:
        U = haar((L, L, 2))
        if b == 0: ms.append(np.mean([action_all(haar((L, L, 2))) for _ in range(3000)])); continue
        U = slice_heatbath(U, b, 30); ss = []
        for i in range(300): U = slice_heatbath(U, b, 1); ss.append(action_all(U))
        ms.append(np.mean(ss))
    out['tempering ladder'] = -np.trapezoid(ms, bs); print("3 tempering ladder", out['tempering ladder'], "%.0fs" % (time.time()-t))
    # 6 Wang Landau energy level sampling for the density of states
    t = time.time(); edges = np.linspace(0, 2*P, 91); nb = len(edges)-1; lg = np.zeros(nb); f = 1.0
    U = haar((L, L, 2)); S = action_all(U); H = np.zeros(nb)
    def binof(s): return min(max(int(np.searchsorted(edges, s)-1), 0), nb-1)
    while f > 1e-3:
        for _ in range(4000):
            x, y, mu = rng.integers(L), rng.integers(L), rng.integers(2)
            old = local_S(U, x, y, mu); keep = U[x, y, mu].copy(); U[x, y, mu] = haar(())
            Sn = S - old + local_S(U, x, y, mu); i0, i1 = binof(S), binof(Sn)
            if np.log(rng.random()) < lg[i0]-lg[i1]: S = Sn; i0 = i1
            else: U[x, y, mu] = keep
            lg[i0] += f; H[i0] += 1
        vis = H > 0
        if vis.sum() > 10 and H[vis].min() > 0.6*H[vis].mean(): f /= 2; H[:] = 0
    vis = lg > 0; lgn = lg - logsumexp(lg[vis])     # normalise: total Haar mass one
    Sc = 0.5*(edges[1:]+edges[:-1])
    out['energy levels (Wang Landau)'] = logsumexp(lgn[vis] - b0*Sc[vis]); print("6 Wang Landau", out['energy levels (Wang Landau)'], "%.0fs" % (time.time()-t))
    # 7 vertical likelihood: g(x) prop p(x) w(S), w = exp(-lg(S)) from the Wang Landau pilot (flat in levels)
    t = time.time(); wl = lambda s: -lg[binof(s)]
    U = haar((L, L, 2)); S = action_all(U); Ss = []
    for i in range(6000):
        U, S, _ = metro_constrained(U, S, None, 18, 0.0, wfun=None) if False else (U, S, 0)
        for _ in range(18):
            x, y, mu = rng.integers(L), rng.integers(L), rng.integers(2)
            old = local_S(U, x, y, mu); keep = U[x, y, mu].copy(); U[x, y, mu] = haar(())
            Sn = S - old + local_S(U, x, y, mu)
            if np.log(rng.random()) < wl(Sn) - wl(S): S = Sn
            else: U[x, y, mu] = keep
        Ss.append(S)
    Ss = np.array(Ss); lw = np.array([wl(s) for s in Ss])
    out['vertical likelihood'] = logsumexp(-b0*Ss - lw) - logsumexp(-lw); print("7 vertical likelihood", out['vertical likelihood'], "%.0fs" % (time.time()-t))
    np.save("/tmp/vlmc_out.npy", out); print({k: round(v-exact, 4) for k, v in out.items()})
