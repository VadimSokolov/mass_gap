# Moderate deviations of the mean plaquette (Rubin and Sethuraman). Two dimensional SU(2): plaquettes are independent with cumulant function
# K(s) = log[ I_1(beta+s)/(beta+s) ] - log[ I_1(beta)/beta ]. Rate I(x) = sup_s (s(u+x) - K(s)); Gaussian rate x^2/(2v); first Cramer term.
# Also: the curvature term of the heat kernel, log(theta/sin theta) = theta^2/6 + ..., reproduces the one loop term of the ladder.
import numpy as np
from scipy.special import ive
from scipy.optimize import brentq
beta = 2.5; u = lambda b: ive(2, b)/ive(1, b); v = 1-3*u(beta)/beta-u(beta)**2
K = lambda s: np.log(ive(1, beta+s)/(beta+s))+s-np.log(ive(1, beta)/beta)       # ive carries exp(-b)
k3 = (lambda h: (u(beta+h)-2*u(beta)+u(beta-h))/h**2)(1e-3)
def rate(x): s = brentq(lambda s: u(beta+s)-(u(beta)+x), -beta+1e-6, 400); return s*(u(beta)+x)-K(s)
print("beta %.1f  mean %.5f  variance %.5f  third cumulant %.5f" % (beta, u(beta), v, k3))
print("  n      alpha   x = 0.3 n^-alpha   n I(x)    ratio to Gaussian   ratio with Cramer term")
for n in (10**3, 10**5, 10**7):
    for a in (0.45, 1/3, 0.2, 0.0):
        x = 0.3*n**(-a); I = rate(x); g = x*x/(2*v); c = g-k3*x**3/(6*v**3); print("%8d  %.2f   %.5f      %9.2f   %.4f              %.4f" % (n, a, x, n*I, I/g, I/c))
for k in (2, 4):
    t = np.array([0.05, 0.1, 0.2]); s = 1/(k/t-(k-1)/6); print("k = %d: curvature matched time s/(t/k) - 1 =" % k, np.round(s/(t/k)-1, 5), " one loop (1-1/k) t/6 =", np.round((1-1/k)*t/6, 5))
