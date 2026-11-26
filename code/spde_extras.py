# Computations for the chapters on stochastic quantisation, filtering and the Langlands volume.
# 1. Spectrum of the Langevin generator of one SU(2) plaquette, L = Laplacian/2 + (beta/2) grad(cos theta) . grad, in the sectors l = 0 (class
#    functions) and l = 1, against the Bakry Emery bound 1 - beta/2.
# 2. The bracket B(t) of the continuous ladder by reprojected Monte Carlo: how many particles are needed.
# 3. Keilson's theorem for the walk on representations of SU(2): passage time from the trivial representation to dimension N.
# 4. One loop coefficient of the ladder per unit of log b.
import numpy as np, sys
from scipy.linalg import eigh
sys.path.insert(0, "code"); from continuous_ladder import logK, bracket, th, haar
def gaps(beta, l, n=1500):
    x = np.linspace(0, np.pi, n+1); mid = 0.5*(x[1:]+x[:-1]); h = x[1]-x[0]; m = np.sin(mid)**2*np.exp(beta*(np.cos(mid)-1))
    K = np.zeros((n+1, n+1)); M = np.zeros(n+1)
    for i in range(n):
        K[i:i+2, i:i+2] += 0.5*m[i]/h*np.array([[1, -1], [-1, 1]])+(l*(l+1)/(2*np.sin(mid[i])**2))*m[i]*h/4*np.ones((2, 2)); M[i:i+2] += m[i]*h/2
    idx = slice(0, n+1) if l == 0 else slice(1, n); ev = eigh(K[idx, idx], np.diag(M[idx]), eigvals_only=True, subset_by_index=[0, 2]); return ev[1] if l == 0 else ev[0]
if __name__ == "__main__":
    out = []
    for beta in (0, 0.5, 1, 2, 3, 4, 6, 8, 12):
        g0, g1 = gaps(beta, 0), gaps(beta, 1); out.append((beta, g0, g1)); print("beta %4.1f  class sector gap %.4f   l = 1 gap %.4f   Bakry Emery %.2f" % (beta, g0, g1, 1-beta/2), flush=True)
    np.save("code/plaquette_gaps.npy", np.array(out))
    rng = np.random.default_rng(8); t = 0.25; L = logK(t, th); p = np.exp(L-L.max())*haar; cdf = np.cumsum(p)/p.sum(); print("bracket B(0.25) exact %.6f" % bracket(t))
    for N in (10**4, 10**5, 10**6):
        est = []
        for r in range(20):
            a = np.interp(rng.uniform(size=N), cdf, th); la = logK(t, a); c = np.cos(a); est.append(t-(4/3)*np.exp(3*t/4)*np.cov(c, la)[0, 1])
        print("  N = %7d  Monte Carlo B = %.6f +- %.6f (sd of one run)" % (N, np.mean(est), np.std(est, ddof=1)))
    for N in (2, 5, 20, 100):
        j = np.arange(1, N); rates = 1-np.cos(j*np.pi/N); P = np.zeros((N-1, N-1))       # states n = 1..N-1, killed at N; n -> n+-1 w.p. (n+-1)/(2n)
        for n in range(1, N):
            if n > 1: P[n-1, n-2] = (n-1)/(2*n)
            if n < N-1: P[n-1, n] = (n+1)/(2*n)
        ev = np.sort(1-np.linalg.eigvals(P).real); mean = np.linalg.solve(np.eye(N-1)-P, np.ones(N-1))[0]
        print("N = %3d  Dirichlet rates match 1 - cos(j pi/N): %s   mean passage time %.4f   (N^2 - 1)/3 = %.4f" % (N, np.allclose(ev, np.sort(rates)), mean, (N*N-1)/3))
    for b in (1.0001, 1.5, 1.93, 2, 3, 4): print("b = %-6g one loop coefficient per log b: %.4f   (true 11/(3 pi^2) = %.4f)" % (b, (b*b-1)/(3*b*b*np.log(b)), 11/(3*np.pi**2)))
