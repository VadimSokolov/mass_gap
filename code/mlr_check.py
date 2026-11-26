"""Monotone likelihood ratio of the SU(2) heat kernel in t: d/dtheta log K_t(theta) increasing in t.
Character series for t >= 2, Poisson resummed image sum for t < 2. Polson and Sokolov, Nov 2026."""
import numpy as np
th = np.linspace(1e-3, np.pi-1e-3, 1500)
def dlogK_series(t):
    n = np.arange(1, 400)[:, None]; T = th[None, :]
    w = n*np.exp(-t*(n**2-1)/4)
    K = (w*np.sin(n*T)).sum(0)/np.sin(th)
    dK = (w*(n*np.cos(n*T)*np.sin(T) - np.sin(n*T)*np.cos(T))).sum(0)/np.sin(th)**2
    return dK/K
def dlogK_images(t):
    m = np.arange(-6, 7)[:, None]; x = th[None, :] + 2*np.pi*m
    e = np.exp(-(x**2 - th[None, :]**2)/t)        # normalised by exp(-theta^2/t)
    S = (x*e).sum(0); dS = ((1-2*x**2/t)*e).sum(0)
    return dS/S - 1/np.tan(th)
print("overlap check at t=2: max diff", np.max(np.abs(dlogK_series(2.0)-dlogK_images(2.0))))
ts = np.geomspace(0.004, 40, 800)
D = np.array([dlogK_images(t) if t < 2 else dlogK_series(t) for t in ts])
dd = np.diff(D, axis=0)
print("MLR holds on 800 x 1500 grid:", np.all(dd > -1e-12), " min increment:", dd.min())
