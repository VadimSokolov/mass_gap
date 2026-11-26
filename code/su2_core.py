# Shared core: SU(2) lattice gauge theory in D dimensions with unit quaternions, checkerboard heat bath, loops, staples.
import numpy as np
rng = np.random.default_rng(1)
def qm(a, b):
    a0,a1,a2,a3 = np.moveaxis(a,-1,0); b0,b1,b2,b3 = np.moveaxis(b,-1,0)
    return np.stack([a0*b0-a1*b1-a2*b2-a3*b3, a0*b1+a1*b0+a2*b3-a3*b2, a0*b2-a1*b3+a2*b0+a3*b1, a0*b3+a1*b2-a2*b1+a3*b0], -1)
cj = lambda a: a*np.array([1.,-1,-1,-1])
sh = lambda A, mu, s: np.roll(A, -s, axis=mu)
def staple(U, mu):
    D = U.shape[-2]; S = 0
    for nu in range(D):
        if nu == mu: continue
        Um, Un = U[..., mu, :], U[..., nu, :]
        S = S + qm(qm(sh(Un, mu, 1), cj(sh(Um, nu, 1))), cj(Un))
        Unm = sh(Un, nu, -1); S = S + qm(qm(cj(sh(Unm, mu, 1)), cj(sh(Um, nu, -1))), Unm)
    return S
def kp(kap):
    out = np.empty_like(kap); todo = np.ones(kap.shape, bool)
    while todo.any():
        k = kap[todo]; r1, r2, r3, r4 = rng.random((4, k.size)); lam2 = -(np.log(r1)+np.cos(2*np.pi*r2)**2*np.log(r3))/(2*k)
        ok = r4**2 <= 1-lam2; idx = np.flatnonzero(todo); out[idx[ok]] = 1-2*lam2[ok]; todo[idx[ok]] = False
    return out
def sweep(U, beta):
    D = U.shape[-2]; par = np.indices(U.shape[:D]).sum(0) % 2
    for mu in range(D):
        for p in (0, 1):
            A = staple(U, mu); k = np.linalg.norm(A, axis=-1); m = par == p; kk = k[m]; V = A[m]/kk[:, None]
            a0 = kp(beta*kk); d = rng.standard_normal((a0.size, 3)); d /= np.linalg.norm(d, axis=1)[:, None]
            W = np.concatenate([a0[:, None], np.sqrt(np.clip(1-a0**2, 0, None))[:, None]*d], 1)
            Uf = U[..., mu, :]; Uf[m] = qm(W, cj(V)); U[..., mu, :] = Uf
    return U
def cold(L, D): U = np.zeros((L,)*D+(D, 4)); U[..., 0] = 1; return U
def line(U, mu, R):
    A = U[..., mu, :]
    for r in range(1, R): A = qm(A, np.roll(U[..., mu, :], -r, axis=mu))
    return A
def loop(U, mu, nu, R, T):
    A, B = line(U, mu, R), line(U, nu, T)
    return qm(qm(qm(A, np.roll(B, -R, axis=mu)), cj(np.roll(A, -T, axis=nu))), cj(B))
