# Particle Migdal Kadanoff: the ladder as a filter recursion on a cloud of SU(2) elements.
# update (bond moving): reweight by f^(k-1) and resample;  predict (decimation): multiply b^2 independent draws.
import numpy as np
rng = np.random.default_rng(4); NP, NC = 200000, 14; n = np.arange(1, NC+1)
def qm(a, b):
    a0,a1,a2,a3 = a.T; b0,b1,b2,b3 = b.T
    return np.stack([a0*b0-a1*b1-a2*b2-a3*b3, a0*b1+a1*b0+a2*b3-a3*b2, a0*b2-a1*b3+a2*b0+a3*b1, a0*b3+a1*b2-a2*b1+a3*b0], 1)
def lift(th):
    d = rng.standard_normal((th.size, 3)); d /= np.linalg.norm(d, axis=1)[:, None]
    return np.column_stack([np.cos(th), np.sin(th)[:, None]*d])
chi = lambda th: np.sin(np.outer(n, th))/np.sin(th)
def heat_angles(beta, size):          # rejection from Haar using the truncated heat kernel series
    u = np.exp(-(2/beta)*(n**2-1)/4); out = np.empty(0)
    fmax = (n*n*u).sum()
    while out.size < size:
        th = np.arccos(1-2*rng.beta(1.5, 1.5, 4*size))
        f = (n*u)@chi(th); out = np.concatenate([out, th[rng.random(th.size)*fmax < f]])
    return out[:size]
def step(th, k=4, b2=4):
    uh = chi(th).mean(1)/n                              # estimated character coefficients
    f = np.clip((n*uh)@chi(th), 1e-12, None); wts = f**(k-1); wts /= wts.sum()
    th = th[rng.choice(th.size, th.size, p=wts)]        # update: bond moving
    g = lift(th)
    for j in range(b2-1): g = qm(g, lift(rng.permutation(th)))   # predict: decimation
    return np.arccos(np.clip(g[:, 0], -1, 1))
def exact(beta, steps, M=4000, NMAX=400):
    x, w = np.polynomial.legendre.leggauss(M); t = (x+1)*np.pi/2; w = w*np.pi/2; haar = (2/np.pi)*np.sin(t)**2*w
    m = np.arange(1, NMAX+1); C = np.sin(np.outer(m, t))/np.sin(t); u = np.exp(-(2/beta)*(m**2-1)/4); out = [u[1]]
    for s in range(steps):
        f = np.clip((m*u)@C, 1e-300, None); g = f**4; c = (C@(g*haar))/m; u = (c/c[0])**4; out.append(u[1])
    return np.array(out)
if __name__ == "__main__":
    res = {}
    for beta in [3.0, 2.0]:
        th = heat_angles(beta, NP); tr = [np.cos(th).mean()]; hist = [th.copy()]
        for s in range(5): th = step(th); tr.append(np.cos(th).mean()); hist.append(th.copy())
        ex = exact(beta, 5); res[beta] = (np.array(tr), ex)
        print("beta", beta, "\n particle u2", np.round(tr, 4), "\n exact    u2", np.round(ex, 4))
        if beta == 3.0: np.save("/tmp/pmk_hist.npy", np.array(hist))
    np.save("/tmp/pmk_tr.npy", np.array([np.concatenate(res[b]) for b in res]))
