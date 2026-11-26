# Adversarial search for a counterexample to the conjecture that convolution on SU(2) preserves the class L (real zeros <= -1 in cos theta).
# Factors with repeated zeros are the case not covered by the theorem: f = prod (1 + x/g_i)^{m_i}, with zeros clustered, spread, or extreme.
import mpmath as mp, numpy as np, sys
sys.path.insert(0, "code"); import real_zeros as R
mp.mp.dps = 160; R.U[:] = R.cheb_u(80)
def build(gam, mult):
    roots = []
    for g, m in zip(gam, mult): roots += [mp.mpf(float(g))]*int(m)
    return R.from_roots(roots)
def check(p):
    while len(p) > 1 and p[-1] == 0: p = p[:-1]
    z = mp.polyroots(p[::-1], maxsteps=3000, extraprec=1500); im = max(abs(mp.im(r))/(1+abs(r)) for r in z); return im < mp.mpf(10)**(-40) and all(mp.re(r) <= -1 for r in z), im, max(mp.re(r) for r in z)
if __name__ == "__main__":
    rng = np.random.default_rng(17); tested = 0; fails = []; closest = -9
    for trial in range(260):
        kind = trial % 5
        def law():
            nf = int(rng.integers(1, 4)); mult = rng.choice([1, 2, 3, 4], nf); scale = [1e-3, 0.05, 1.0, 30.0, 1e3][int(rng.integers(0, 5))]
            gam = 1+rng.exponential(scale, nf) if kind else 1+np.array([1e-4, 1e4, 3.0])[:nf]; return gam, mult
        (ga, ma), (gb, mb) = law(), law()
        if sum(ma)+sum(mb) > 22: continue
        ok, im, top = check(R.conv(build(ga, ma), build(gb, mb))); tested += 1; closest = max(closest, float(top))
        if not ok: fails.append((ga, ma, gb, mb, mp.nstr(im, 3)))
    print("pairs with repeated zeros tested: %d;  outside the class: %d;  largest zero seen: %.8f" % (tested, len(fails), closest))
    for f in fails[:5]: print("  FAIL", f)
    # the first rung from a heat kernel, with the product truncated to M factors
    from scipy.optimize import brentq
    for t, M in ((2.0, 5), (1.0, 6)):
        N = 400; C = np.cosh((2*np.arange(1, N+1)-1)*t/4); f = lambda g: np.sum(1/(C-g)); gam = [brentq(f, C[m]*(1+1e-13)+1e-13, C[m+1]*(1-1e-13)-1e-13) for m in range(M)]
        rho = build(gam, [4]*M); c2 = R.conv(rho, rho); ok, im, top = check(R.conv(c2, c2)); print("heat kernel t = %.1f, first %d factors, power 4 then fourfold convolution: in class %s, largest zero %s (largest zero before: %.5f)" % (t, M, ok, mp.nstr(top, 8), -gam[0]))
