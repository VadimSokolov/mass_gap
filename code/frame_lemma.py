"""Frame lemma (Chapter on theta kernels): random checks of
   rank Q >= (tr Q)^2 / tr Q^2  and
   codim sum_i B_i U_i <= (eta + 4 alpha + 4 sqrt(e)) D,
   and of the layer cake bound for a mean."""
import numpy as np
rng = np.random.default_rng(7)

def trial(D, nmaps, dv, drop, noise):
    Bs, Us = [], []
    for i in range(nmaps):
        G = rng.standard_normal((D, dv)) + 1j*rng.standard_normal((D, dv))
        Qm, _ = np.linalg.qr(G)                       # isometry V_i -> H
        Qm = Qm @ (np.eye(dv) + noise*rng.standard_normal((dv, dv)))
        Bs.append(Qm)
        Us.append(np.eye(dv)[:, :dv-drop])            # U_i: drop some directions
    c = np.ones(nmaps)
    S = nmaps*dv
    Q = sum(ci*B @ B.conj().T for ci, B in zip(c, Bs))
    e = sum(np.linalg.norm(ci*B.conj().T @ B - np.eye(dv))**2 for ci, B in zip(c, Bs))/S
    trQ = np.trace(Q).real
    eta = max(np.trace(Q @ Q).real*D/trQ**2 - 1, 0)
    alpha = nmaps*drop/S
    W = np.hstack([B @ U for B, U in zip(Bs, Us)])
    codim = D - np.linalg.matrix_rank(W, tol=1e-9)
    rankQ = np.linalg.matrix_rank(Q, tol=1e-9)
    rho = (alpha + np.sqrt(alpha*e))/(1 - np.sqrt(e)) if e < 1 else np.inf
    sharp = 0.25*(np.sqrt(eta) + np.sqrt(eta + 4*rho))**2*D     # sharp form of the lemma
    return codim, (eta + 4*alpha + 4*np.sqrt(e))*D, rankQ, trQ**2/np.trace(Q@Q).real, e, sharp

bad = 0; n = 0; tight = 0; tight2 = 0; worse = 0
for D in (12, 20, 40):
    for nmaps in (1, 2, 3, 5, 8):
        for dv in (2, 4, 6):
            for drop in (0, 1):
                for noise in (0.0, 0.05, 0.15):
                    cd, bound, rk, lower, e, sharp = trial(D, nmaps, dv, drop, noise)
                    n += 1
                    if e <= 0.25 and cd > bound + 1e-9: bad += 1
                    if rk < lower - 1e-9: bad += 1
                    if e <= 0.25 and bound < D: tight += 1
                    if e < 1 and cd > sharp + 1e-9: bad += 1
                    if e < 1 and sharp < D: tight2 += 1
                    if e <= 0.25 and sharp > bound + 1e-9: worse += 1
print("trials", n, "violations", bad, "trials with a bound below D", tight,
      "with the sharp bound below D", tight2, "sharp worse than plain", worse)

# a tight frame: an orthogonal decomposition of H, where codim = alpha D exactly
D = 24; dv = 6; Bs = [np.eye(D)[:, i*dv:(i+1)*dv] for i in range(4)]
W = np.hstack([B[:, :dv-1] for B in Bs])
print("orthogonal decomposition: codim", D - np.linalg.matrix_rank(W), "alpha D", 4*1/24*D)

# layer cake: mean >= a - eps0*H - 4u log(1+H/u) when N(t)/D <= eps0 + 4u/(u+t)
viol = 0
for _ in range(2000):
    D = 200; a = 0.0; H = rng.uniform(1, 10); u = rng.uniform(0.01, 1); eps0 = rng.uniform(0, 0.2)
    d = np.clip(a - rng.exponential(u, D)*rng.integers(0, 2, D), a - H, None) + rng.uniform(0, 1, D)*0
    ts = np.linspace(0, H, 400)
    N = np.array([(d < a - t).sum() for t in ts])/D
    if np.all(N <= eps0 + 4*u/(u+ts)):
        if d.mean() < a - eps0*H - 4*u*np.log(1+H/u) - 1e-12: viol += 1
print("layer cake violations", viol)
