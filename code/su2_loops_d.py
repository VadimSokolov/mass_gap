# Laws of R x R Wilson loop holonomies, R = 1, 2, 4, for SU(2) in D dimensions. Usage: beta D L
import numpy as np, sys
rng = np.random.default_rng(3); D = int(sys.argv[2]); L = int(sys.argv[3])
def qm(a, b):
    a0,a1,a2,a3 = np.moveaxis(a,-1,0); b0,b1,b2,b3 = np.moveaxis(b,-1,0)
    return np.stack([a0*b0-a1*b1-a2*b2-a3*b3, a0*b1+a1*b0+a2*b3-a3*b2, a0*b2-a1*b3+a2*b0+a3*b1, a0*b3+a1*b2-a2*b1+a3*b0], -1)
cj = lambda a: a*np.array([1.,-1,-1,-1])
sh = lambda A, mu, s: np.roll(A, -s, axis=mu)          # A(x + s mu)
def staple(U, mu):
    S = 0
    for nu in range(D):
        if nu == mu: continue
        Um, Un = U[..., mu, :], U[..., nu, :]
        S = S + qm(qm(sh(Un, mu, 1), cj(sh(Um, nu, 1))), cj(Un))
        Unm = sh(Un, nu, -1)
        S = S + qm(qm(cj(sh(Unm, mu, 1)), cj(sh(Um, nu, -1))), Unm)
    return S
def kp(kap):
    """sample a0 with density sqrt(1-a0^2) exp(kap a0), Kennedy Pendleton, vectorised"""
    out = np.empty_like(kap); todo = np.ones(kap.shape, bool)
    while todo.any():
        k = kap[todo]; r1, r2, r3, r4 = rng.random((4, k.size))
        lam2 = -(np.log(r1) + np.cos(2*np.pi*r2)**2*np.log(r3))/(2*k)
        ok = r4**2 <= 1-lam2
        idx = np.flatnonzero(todo); out[idx[ok]] = 1-2*lam2[ok]; todo[idx[ok]] = False
    return out
def sweep(U, beta, par):
    for mu in range(D):
        A = staple(U, mu); k = np.linalg.norm(A, axis=-1)
        for p in (0, 1):
            m = par == p; kk = k[m]; V = A[m]/kk[:, None]
            a0 = kp(beta*kk); d = rng.standard_normal((a0.size, 3)); d /= np.linalg.norm(d, axis=1)[:, None]
            W = np.concatenate([a0[:, None], np.sqrt(np.clip(1-a0**2, 0, None))[:, None]*d], 1)
            Unew = qm(W, cj(V)); Ufull = U[..., mu, :]; Ufull[m] = Unew; U[..., mu, :] = Ufull
            A = staple(U, mu); k = np.linalg.norm(A, axis=-1)
    return U
def line(U, mu, R):
    A = U[..., mu, :]
    for r in range(1, R): A = qm(A, np.roll(U[..., mu, :], -r, axis=mu))
    return A
def loop(U, mu, nu, R):
    A, B = line(U, mu, R), line(U, nu, R)
    return qm(qm(qm(A, np.roll(B, -R, axis=mu)), cj(np.roll(A, -R, axis=nu))), cj(B))
beta = float(sys.argv[1]); g = np.indices((L,)*D).sum(0) % 2
U = np.zeros((L,)*D+(D, 4)); U[..., 0] = 1
for i in range(100): U = sweep(U, beta, g)
bins = np.linspace(0, np.pi, 61); H = {1: np.zeros(60), 2: np.zeros(60), 4: np.zeros(60)}; W = {1: [], 2: [], 4: []}
for i in range(80):
    for j in range(4): U = sweep(U, beta, g)
    for mu in range(D):
        for nu in range(mu+1, D):
            for R in (1, 2, 4):
                a0 = loop(U, mu, nu, R)[..., 0]; H[R] += np.histogram(np.arccos(np.clip(a0, -1, 1)), bins)[0]; W[R].append(a0.mean())
np.save("/tmp/loopsd_%d_%s.npy" % (D, sys.argv[1]), np.array([H[1], H[2], H[4]]))
print("D=%d beta %.1f  W11 %.5f W22 %.5f W44 %.5f" % (D, beta, np.mean(W[1]), np.mean(W[2]), np.mean(W[4])), flush=True)
