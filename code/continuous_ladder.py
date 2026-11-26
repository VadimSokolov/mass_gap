# The ladder in continuous scale, as a filtering equation. With b = 1 + dl the bond moving step is a Bayes update with log likelihood
# (d-2) dl log rho and decimation is a prediction step of time 2 dl t. For the heat kernel projection (u_2 = exp(-3t/4)):
#     dt/dl = (4-d) t + (d-2) [ t - (4/3) exp(3t/4) Cov_{K_t}(cos theta, log K_t) ].
# Part 1 evaluates the bracket and compares with the series t^2/6 + t^3/36 + t^4/432. Part 2 is a weighted particle system of
# Kurtz and Xiong type for the same flow in d = 4: Brownian particles on SU(2) with weights driven by the log density.
import numpy as np
x, w = np.polynomial.legendre.leggauss(4000); th = (x+1)*np.pi/2; haar = (2/np.pi)*np.sin(th)**2*w*np.pi/2
def logK(t, a):
    if t >= 1.0:
        n = np.arange(1, 200)[:, None]; return np.log(np.clip((n*np.exp(-(n*n-1)*t/4)*np.sin(n*a)/np.sin(a)).sum(0), 1e-300, None))
    j = np.arange(-4, 5)[:, None]; aj = a+2*np.pi*j; m = (-(aj**2)/t).max(0); return m+np.log(np.clip((aj*np.exp(-(aj**2)/t-m)).sum(0)/np.sin(a), 1e-300, None))
def bracket(t):
    L = logK(t, th); p = np.exp(L-L.max())*haar; p /= p.sum(); c = np.cos(th); return t-(4/3)*np.exp(3*t/4)*((p*c*L).sum()-(p*c).sum()*(p*L).sum())
if __name__ == "__main__":
    print("   t     bracket      series t^2/6+t^3/36+t^4/432   d beta/dl in d=4 (beta = 2/t)")
    for t in (0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4):
        B = bracket(t); print("%5.2f  %.8f   %.8f                  %.5f" % (t, B, t*t/6+t**3/36+t**4/432, -4*B/t**2))
    rng = np.random.default_rng(4); N = 20000; dl = 0.01; t = 0.25  # hybrid system, no reprojection
    q = np.zeros((N, 4)); q[:, 0] = 1; s = 0.0
    while s < t/2-1e-12:                      # start from K_t: Brownian motion with generator Laplacian/4 run for time t
        xi = rng.standard_normal((N, 4))*np.sqrt(0.001); xi -= (xi*q).sum(1, keepdims=True)*q; q += xi; q /= np.linalg.norm(q, axis=1, keepdims=True); s += 0.001
    tt = -(4/3)*np.log(q[:, 0].mean()); ode = tt; out = []
    for step in range(301):
        if step % 50 == 0: out.append((step*dl, tt, ode)); print("l = %.2f   particle t = %.4f   ode t = %.4f   beta = %.3f" % (step*dl, tt, ode, 2/tt), flush=True)
        a = np.arccos(np.clip(q[:, 0], -1, 1)); a = np.clip(a, 1e-6, np.pi-1e-6); lw = 2*dl*logK(tt, a)           # update: weight by rho^{(d-2) dl}
        wgt = np.exp(lw-lw.max()); q = q[rng.choice(N, N, p=wgt/wgt.sum())]                                  # resample
        nsub = max(1, int(np.ceil(dl*tt/0.002))); h = 2*dl*tt/2/nsub                                         # predict: heat flow for time 2 dl t
        for _ in range(nsub):
            xi = rng.standard_normal((N, 4))*np.sqrt(h); xi -= (xi*q).sum(1, keepdims=True)*q; q += xi; q /= np.linalg.norm(q, axis=1, keepdims=True)
        tt = -(4/3)*np.log(q[:, 0].mean()); ode += dl*2*bracket(ode)
    np.save("code/continuous_ladder.npy", np.array(out))
