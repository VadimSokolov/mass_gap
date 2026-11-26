"""Numerical checks for the Thorin / Yang Mills deck and note.
Polson and Sokolov, Nov 2026."""
import numpy as np, mpmath as mp
from scipy.integrate import quad
mp.mp.dps = 30

# 1. The 2 sin^2 lemma: e^{-wx}+wx-1 = int_0^inf log(1+u/z) tau_x(z) dz, u=w^2/2
def lemma_rhs(x, w):
    u = w*w/2
    # substitute z = 2 y^2: integrand (4x/pi) sin^2(xy) log(1+w^2/(4y^2)) dy
    f = lambda y: (4*x/mp.pi)*mp.sin(x*y)**2*mp.log(1+w**2/(4*y**2))
    return mp.quadosc(f, [0, mp.inf], omega=2*x)
for x, w in [(0.7, 1.3), (2.0, 0.5), (3.1, 2.2)]:
    lhs = mp.e**(-w*x)+w*x-1
    print("lemma", x, w, mp.nstr(lhs, 12), mp.nstr(lemma_rhs(x, w), 12))

# 2. Toy glueball spectral measure: atoms + continuum from 2*Delta
Delta = 1.0
atoms = [(1.0, 1.0), (1.5, 0.6), (1.9, 0.3)]           # (E, weight)
def cont_density(E):                                    # two particle continuum
    return np.where(E > 2*Delta, 0.4*np.sqrt(np.clip(E-2*Delta, 0, None))*np.exp(-0.3*(E-2*Delta)), 0.0)
def C(tau):  # complex tau allowed
    s = sum(c*np.exp(-tau*E) for E, c in atoms)
    s += quad(lambda E: (cont_density(E)*np.exp(-tau*E)).real, 2*Delta, 80)[0] if np.isreal(tau) else \
         quad(lambda E: cont_density(E)*np.exp(-tau.real*E)*np.cos(tau.imag*E), 2*Delta, 80, limit=400)[0]
    return s
def T_direct(alpha, z):
    y = np.sqrt(z/2); t = np.sqrt(2*z)
    s = sum(c*np.exp(-alpha*E)*2*np.sin(E*y)**2 for E, c in atoms)
    s += quad(lambda E: cont_density(E)*np.exp(-alpha*E)*2*np.sin(E*y)**2, 2*Delta, 80, limit=400)[0]
    return s/(np.pi*t)
def T_realtime(alpha, z):
    t = np.sqrt(2*z)
    return (C(alpha) - C(complex(alpha, t)))/(np.pi*t)
for a, z in [(0.5, 0.8), (1.0, 3.0), (2.0, 10.0)]:
    print("T formula", a, z, T_direct(a, z), T_realtime(a, z))

# 3. Gap from Thorin decay rate
for z in [0.8, 3.0]:
    for a in [5, 10, 20, 40]:
        print("rate z=%g alpha=%g  -(1/alpha) log T = %.4f" % (z, a, -np.log(T_direct(a, z))/a))

# 4. Riemann kernel from the SU(2) sphere partition function
def Phi(u, N=60):
    return sum((2*mp.pi**2*n**4*mp.e**(4.5*u) - 3*mp.pi*n**2*mp.e**(2.5*u))*mp.e**(-mp.pi*n**2*mp.e**(2*u)) for n in range(1, N))
def F(A, N=60):   # F(A) = sum n^2 e^{-A n^2/8} = e^{-A/8} Z_0(A)
    return sum(n**2*mp.e**(-A*n**2/8) for n in range(1, N))
def Fp(A, N=60):
    return sum(-n**4/8*mp.e**(-A*n**2/8) for n in range(1, N))
for u in [0.0, 0.3, -0.2]:
    A = 8*mp.pi*mp.e**(2*u)
    print("Phi", u, mp.nstr(Phi(u), 15), mp.nstr(2*mp.pi**2*mp.e**(4.5*u)*(-8*Fp(A)) - 3*mp.pi*mp.e**(2.5*u)*F(A), 15))

# 5. 2D SU(2) torus: gap 3/8 from the Thorin rate of Z_1 minus vacuum
def T2d(A0, z, g=1, N=80):
    y = np.sqrt(z/2); t = np.sqrt(2*z)
    s = sum(n**(2-2*g)*np.exp(-A0*(n*n-1)/8)*2*np.sin((n*n-1)/8*y)**2 for n in range(2, N))
    return s/(np.pi*t)
for A0 in [10, 40, 80]:
    print("2D torus rate A0=%g: %.5f (gap 3/8 = 0.375)" % (A0, -np.log(T2d(A0, 5.0))/A0))
