# (a) One plaquette in the Hamiltonian formulation: H = lam * Casimir + (1 - cos theta)/lam on class functions of SU(2),
#     in the character basis chi_1..chi_N, where the Casimir is diag((n^2-1)/4) and cos theta is tridiagonal with entries 1/2.
#     Gap, the number of levels needed, and the count E[1/(1-cos theta)] - 2 of the ground state.
# (b) A band limited model for GBC: pooled characters of n = 256 Wilson plaquettes to band N, a cubic polynomial readout fitted by ridge
#     regression (a convex problem), against the exact posterior median of beta.
import numpy as np, sys
def ham(lam, N): n = np.arange(1, N+1); return np.diag(lam*(n*n-1)/4+1/lam)-np.diag(np.full(N-1, 0.5/lam), 1)-np.diag(np.full(N-1, 0.5/lam), -1)
def gap(lam, N): e = np.linalg.eigvalsh(ham(lam, N)); return e[1]-e[0]
def levels(lam, tol=1e-8):
    ref = gap(lam, 600); N = 2
    while abs(gap(lam, N)-ref) > tol: N += 1
    return N
def count(lam, N=600):
    e, v = np.linalg.eigh(ham(lam, N)); th = (np.arange(200000)+0.5)*np.pi/200000; psi = sum(v[k, 0]*np.sin((k+1)*th)/np.sin(th) for k in range(min(N, 200))); w = psi**2*np.sin(th)**2
    return (w/(1-np.cos(th))).sum()/w.sum()-2
if __name__ == "__main__":
    print("%8s %10s %10s %8s %10s %12s" % ("lambda", "gap", "3 lam/4", "levels", "count", "4 sqrt2/lam"))
    for lam in (8.0, 2.0, 1.0, 0.5, 0.2, 0.1, 0.05, 0.02):
        print("%8.2f %10.5f %10.5f %8d %10.3f %12.3f" % (lam, gap(lam, 600), 0.75*lam, levels(lam), count(lam), 4*np.sqrt(2)/lam))
    print("weak coupling limit of the gap: sqrt(2) = %.5f" % np.sqrt(2))
    sys.path.insert(0, "code"); sys.path.insert(0, "."); import gbc_exact as G
    G.rng = np.random.default_rng(21); n = G.n
    def embed(a, B): th = np.clip(np.arccos(np.clip(a, -1, 1)), 1e-9, np.pi-1e-9); return np.array([(np.sin(k*th)/np.sin(th)).mean()/k for k in range(2, B+2)])
    def sim(M, B=8): b = G.rng.uniform(0.5, 6, M); return b, np.array([embed(G.draw(x, n), B) for x in b])
    def poly(S):
        cols = [np.ones(len(S))]+[S[:, i] for i in range(S.shape[1])]+[S[:, i]*S[:, j] for i in range(S.shape[1]) for j in range(i, S.shape[1])]+[S[:, 0]**3]; return np.column_stack(cols)
    b, S = sim(4000); bt, St = sim(500); ex = np.array([G.exact_q(s[0])[2] for s in St]); mu, sd = S.mean(0), S.std(0)
    for B in (1, 2, 4, 8):
        X = poly((S[:, :B]-mu[:B])/sd[:B]); Xt = poly((St[:, :B]-mu[:B])/sd[:B]); w = np.linalg.solve(X.T@X+1e-3*np.eye(X.shape[1]), X.T@b)
        print("band %d: %3d features, mean abs error of the readout against the exact posterior median %.4f" % (B, X.shape[1], np.abs(Xt@w-ex).mean()))
