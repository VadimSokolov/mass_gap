# Vectorised checkerboard heat bath for 4D SU(2) Wilson lattice gauge theory; measures the mean plaquette.
import numpy as np, sys
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
def qm(a, b):
    a0,a1,a2,a3 = np.moveaxis(a,-1,0); b0,b1,b2,b3 = np.moveaxis(b,-1,0)
    return np.stack([a0*b0-a1*b1-a2*b2-a3*b3, a0*b1+a1*b0+a2*b3-a3*b2, a0*b2-a1*b3+a2*b0+a3*b1, a0*b3+a1*b2-a2*b1+a3*b0], -1)
cj = lambda a: a*np.array([1.,-1,-1,-1])
sh = lambda A, mu, s: np.roll(A, -s, axis=mu)          # A(x + s mu)
def staple(U, mu):
    S = 0
    for nu in range(4):
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
    for mu in range(4):
        A = staple(U, mu); k = np.linalg.norm(A, axis=-1)
        for p in (0, 1):
            m = par == p; kk = k[m]; V = A[m]/kk[:, None]
            a0 = kp(beta*kk); d = rng.standard_normal((a0.size, 3)); d /= np.linalg.norm(d, axis=1)[:, None]
            W = np.concatenate([a0[:, None], np.sqrt(np.clip(1-a0**2, 0, None))[:, None]*d], 1)
            Unew = qm(W, cj(V)); Ufull = U[..., mu, :]; Ufull[m] = Unew; U[..., mu, :] = Ufull
            A = staple(U, mu); k = np.linalg.norm(A, axis=-1)
    return U
def plaq(U):
    tot = 0
    for mu in range(4):
        for nu in range(mu+1, 4):
            Um, Un = U[..., mu, :], U[..., nu, :]
            tot += qm(qm(qm(Um, sh(Un, mu, 1)), cj(sh(Um, nu, 1))), cj(Un))[..., 0].mean()
    return tot/6
L = 8; beta = float(sys.argv[1])
g = np.indices((L,)*4).sum(0) % 2
U = np.zeros((L,)*4+(4, 4)); U[..., 0] = 1
for i in range(60): U = sweep(U, beta, g)
P = []
for i in range(600):
    U = sweep(U, beta, g); P.append(plaq(U))
P = np.array(P); bl = P.reshape(30, 20).mean(1)
print("%.3f %.6f %.6f" % (beta, P.mean(), bl.std()/np.sqrt(30)))
np.save("/tmp/plaq_hist_%s.npy" % sys.argv[1], P)
