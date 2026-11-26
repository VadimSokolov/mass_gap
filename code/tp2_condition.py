# Condition C for a class function law f on SU(2): theta -> (Delta f)/f is nondecreasing, Delta = d^2 + 2 cot d the radial Laplacian.
# Equivalently h''/h nondecreasing for h = sin(theta) f. Under C, convolution with f is totally positive of order two.
# We test C for heat kernels, their powers, and the laws along the ladder, and test the TP2 property of the kernel directly.
import numpy as np
from scipy.special import ive
th = np.linspace(0.02, np.pi-0.02, 1500); n = np.arange(1, 300)
SIN = np.sin(np.outer(n, th))
def h_and_derivs(u):                      # h = sin(theta) f = sum n u_n sin(n theta)
    return (n*u)@SIN, (n**3*u)@SIN        # h, -h''
def cond(u):
    h, mh2 = h_and_derivs(u); q = -mh2/h; return np.gradient(q, th)/np.sin(th)
heat = lambda t: np.exp(-(n**2-1)*t/4)
def power(u, k, M=4000):
    x, w = np.polynomial.legendre.leggauss(M); t = (x+1)*np.pi/2; w = w*np.pi/2; C = np.sin(np.outer(n, t))/np.sin(t)
    f = np.clip((n*u)@C, 1e-300, None); c = (C@(f**k*(2/np.pi)*np.sin(t)**2*w))/n; return c/c[0]
print("heat kernels: min of d/dtheta(Delta K/K)/sin over theta")
for t in [0.3, 0.6, 1, 2, 4]: print("  t=%.1f  K: %.4f   K^2: %.4f   K^4: %.4f" % (t, cond(heat(t)).min(), cond(power(heat(t), 2)).min(), cond(power(heat(t), 4)).min()))
print("ladder from Wilson beta = 2.5 (after decimation, and after bond moving)")
u = np.clip(ive(n, 2.5)/ive(1, 2.5), 1e-300, None)
for m in range(6):
    ub = power(u, 4); print("  rung %d: law %.4f   bond moved %.4f" % (m, cond(u).min(), cond(ub).min())); u = np.clip(ub, 1e-300, None)**4
print("ladder from heat kernel beta = 6")
u = heat(2/6.0)
for m in range(6):
    ub = power(u, 4); print("  rung %d: law %.4f   bond moved %.4f" % (m, cond(u).min(), cond(ub).min())); u = np.clip(ub, 1e-300, None)**4
# direct test of total positivity: Psi(a,b) = (int_a^b h)(h'(a)+h'(b)) - h(b)^2 + h(a)^2 >= 0 for 0 < a < b < pi
COS = np.cos(np.outer(n, th))
def psi_min(u):
    h = (n*u)@SIN; hp = (n*n*u)@COS; H = -((u)@COS)            # antiderivative of h: -sum u_n cos(n theta)
    P = (H[None, :]-H[:, None])*(hp[:, None]+hp[None, :])-h[None, :]**2+h[:, None]**2          # P[a, b]
    iu = np.triu_indices(th.size, 5); scale = (h[None, :]**2+h[:, None]**2)[iu]; r = P[iu]/scale
    j = r.argmin(); return r.min(), th[iu[0][j]], th[iu[1][j]], (r < -1e-9).mean()
print("direct TP2 test: min of Psi/(h(a)^2+h(b)^2), where attained, share of pairs violating")
for name, u in [("K_2", heat(2.0)), ("K_2^4", power(heat(2.0), 4)), ("K_1^4", power(heat(1.0), 4)), ("K_0.6^4", power(heat(0.6), 4)), ("Wilson 2.5", np.clip(ive(n, 2.5)/ive(1, 2.5), 1e-300, None)), ("Wilson 2.5 ^4", power(np.clip(ive(n, 2.5)/ive(1, 2.5), 1e-300, None), 4)), ("Wilson 1.4", np.clip(ive(n, 1.4)/ive(1, 1.4), 1e-300, None))]:
    print("  %-14s %.3e at (%.2f, %.2f)  violating %.4f" % ((name,)+psi_min(u)))
