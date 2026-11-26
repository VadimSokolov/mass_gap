"""Migdal Kadanoff flow for 4D SU(2) lattice gauge theory in character (Thorin level) variables.
Plaquette weight f(theta) = sum_n n c_n chi_n(theta), chi_n = sin(n theta)/sin(theta), Haar (2/pi) sin^2.
Step with scale b=2 in d=4: bond moving f -> f^(b^(d-2)) = f^4, decimation c_n -> c_n^(b^2) = c_n^4.
Effective beta: solve I_2(beta)/I_1(beta) = u_2. Polson and Sokolov, Nov 2026."""
import numpy as np
from scipy.special import ive
from scipy.optimize import brentq
M = 6000; NMAX = 1500
x, w = np.polynomial.legendre.leggauss(M)
th = (x+1)*np.pi/2; w = w*np.pi/2; haar = (2/np.pi)*np.sin(th)**2*w
n = np.arange(1, NMAX+1)
CHI = np.sin(np.outer(n, th))/np.sin(th)            # NMAX x M
def coeffs(logf):                                    # c_n / c_1 from log f (stabilised)
    f = np.exp(logf-logf.max())
    c = (CHI@(f*haar))/n
    return c/c[0]
def logf_from(u):                                    # f = sum n u_n chi_n, normalised c_1 = 1
    f = (n*u)@CHI
    return np.log(np.clip(f, 1e-300, None))
def beta_eff(u2):
    return brentq(lambda b: ive(2, b)/ive(1, b)-u2, 1e-8, 1e5)
def step(u, k=4, b2=4):
    lf = logf_from(u)
    u = coeffs(k*lf)                                 # bond moving
    return np.clip(u, 0, None)**b2                   # decimation (exact, Thorin area shift)
def flow(beta0, steps=40, stop=0.02):
    u = ive(n, beta0)/ive(1, beta0)
    bs = [beta0]
    for s in range(steps):
        u = step(u)
        if u[1] < stop: bs.append(beta_eff(max(u[1], 1e-12))); break
        bs.append(beta_eff(u[1]))
    return np.array(bs)
if __name__ == "__main__":
    # check coefficients reproduce Bessel ratios
    b = 3.0; th0 = b*np.cos(th); u = coeffs(th0)
    print("check u_2, u_3:", u[1], ive(2, b)/ive(1, b), u[2], ive(3, b)/ive(1, b))
    for b0 in [0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0]:
        f = flow(b0); print("beta0=%5.1f steps to strong coupling=%2d  flow:" % (b0, len(f)-1), np.round(f[:8], 3))
