"""Open problem 1 checks. Polson and Sokolov, Nov 2026."""
import numpy as np
from scipy.special import ive
from scipy.optimize import brentq
b = np.logspace(-3, 4, 4000)
u2 = ive(2, b)/ive(1, b)                       # I_2/I_1, exponentially scaled ratio is exact
up = b/(1.5+np.sqrt(b**2+2.25))                # claimed Amos upper bound
lo = b/(2+np.sqrt(b**2+4))                     # Amos lower bound
print("upper bound holds:", np.all(u2 <= up*(1+1e-14)), " lower holds:", np.all(u2 >= lo*(1-1e-14)))
print("check identity -log(up) = arcsinh(3/(2b)):", np.max(np.abs(-np.log(up)-np.arcsinh(1.5/b))))
G = -b*np.log(u2)
print("G increasing:", np.all(np.diff(G) > 0), " G(inf) ~", G[-1], "(3/2)")
print("u_n <= u_2 for n=3..10 at all b:", all(np.all(ive(n, b)/ive(1, b) <= u2) for n in range(3, 11)))
# continuum ratio Delta(a)/Delta_cont = G(beta)/(3/2), x = 3/(2 beta)
for beta in [0.5, 1, 2, 5, 20, 100]:
    x = 1.5/beta; r = (-beta*np.log(ive(2, beta)/ive(1, beta)))/1.5
    print("beta=%6.1f x=%.3f  Delta(a)/Delta_cont=%.5f  bound asinh(x)/x=%.5f" % (beta, x, r, np.arcsinh(x)/x))
# 4D strong coupling leading order glueball mass for SU(2): m = -4 log u
for beta in [0.1, 0.5, 1.0, 1.5]:
    print("4D SU(2) beta=%.1f  u=%.4f  m_leading=-4 log u=%.3f" % (beta, ive(2, beta)/ive(1, beta), -4*np.log(ive(2, beta)/ive(1, beta))))
# Open problem 2: zero of psi between atoms -> negative Thorin atom.
# psi_alpha(s) = int mu_alpha(dE)/(E+s), mu_alpha(dE) = e^{-alpha E} mu(dE), for the full toy spectrum of the book (code/figures.py, checks.py):
# three atoms plus the continuum 0.4 sqrt(E-2) exp(-0.3(E-2)) on E > 2 (Delta = 1). The first line (alpha = 0, atoms only) is what the book printed before.
from scipy.integrate import quad
atoms = [(1.0, 1.0), (1.5, 0.6), (1.9, 0.3)]
rho = lambda E: 0.4*np.sqrt(E-2)*np.exp(-0.3*(E-2))
def psi(s, alpha, continuum):
    v = sum(w*np.exp(-alpha*E)/(E+s) for E, w in atoms)
    if continuum: v += quad(lambda E: rho(E)*np.exp(-alpha*E)/(E+s), 2, 200, limit=400)[0]
    return v
for alpha, cont in [(0, False), (0, True), (0.5, True)]:
    zz = [-brentq(lambda s: psi(s, alpha, cont), -E2+1e-9, -E1-1e-9, xtol=1e-14) for (E1, _), (E2, _) in zip(atoms[:-1], atoms[1:])]
    print("alpha=%.1f continuum=%-5s: negative Thorin atoms at %.5f and %.5f" % (alpha, cont, zz[0], zz[1]))
