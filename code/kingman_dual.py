# Kingman's coalescent inside the SU(2) heat kernel. In y = cos^2(theta/2) the radial part of Brownian motion on SU(2) is a neutral
# Wright Fisher diffusion with mutation parameters (3/2, 3/2); its dual is the pure death process A with rates m(m+2)/2 (the number of
# lineages of Kingman's coalescent with mutation theta = 3). Checks:
#  1. The latent count of the heat kernel K_t has the law of A(t/2) started from infinity (Griffiths; Tavare).
#  2. Heat flow acts on the densities f_a by the death process: f_a * K_s = sum_j P_a(A(s/2) = j) f_j.
import numpy as np
from scipy.linalg import expm
from scipy.special import gammaln
from scipy.optimize import brentq
MM = 60; j = np.arange(MM+1); logZ = gammaln(2*j+3)-2*gammaln(j+2)-np.log(j+2)-j*np.log(2)
def u(n, J): return np.exp(gammaln(J)+gammaln(J+2)-gammaln(J-n+1)-gammaln(J+n+1)) if n <= J else 0.0
def count_law(t):                                   # from the product formula, as in latent_mixture.py
    N = int(60/t)+60; C = np.cosh((2*np.arange(1, N+1)-1)*t/4); f = lambda g: np.sum(1/(C-g)); nu = np.array([1.0])
    for m in range(N-40): p = 1/brentq(f, C[m]*(1+1e-13)+1e-13, C[m+1]*(1-1e-13)-1e-13); nu = np.convolve(nu, [1-p, p])[:MM+1]
    nu = np.pad(nu, (0, MM+1-nu.size)); p = nu*np.exp(logZ); return p/p.sum()
def tavare(t, theta=3.0, K=400):                    # P(A(tau) = m), tau = t/2, entrance law from infinity
    tau = t/2; out = []
    for m in range(MM+1):
        k = np.arange(max(m, 1), K); lg = gammaln(m+theta+k-1)-gammaln(m+theta)-gammaln(m+1)-gammaln(k-m+1) if m > 0 else gammaln(theta+k-1)-gammaln(theta)-gammaln(k+1)
        term = np.exp(-k*(k+theta-1)*tau/2+lg)*(2*k+theta-1)*(-1.0)**(k-m); out.append(term.sum() if m > 0 else 1-np.sum(np.exp(-k*(k+theta-1)*tau/2+lg)*(2*k+theta-1)*(-1.0)**(k-1)))
    return np.array(out)
if __name__ == "__main__":
    G = np.zeros((MM+1, MM+1))
    for m in range(1, MM+1): G[m, m] = -m*(m+2)/2; G[m, m-1] = m*(m+2)/2
    for t in (2.0, 1.0, 0.5):
        p = count_law(t); q = tavare(t); print("t %.1f  latent count law vs coalescent entrance law: max abs difference %.2e   means %.4f %.4f" % (t, np.abs(p-q).max(), (p*j).sum(), (q*j).sum()))
    for a, s in ((5, 0.4), (12, 0.15), (30, 0.05)):
        P = expm(G*s/2)[a]; lhs = np.array([u(n, a+1)*np.exp(-(n*n-1)*s/4) for n in range(1, a+2)]); rhs = np.array([sum(P[i]*u(n, i+1) for i in range(a+1)) for n in range(1, a+2)])
        print("f_%d * K_%.2f: characters from the death process agree to %.1e;  mean count %.3f" % (a, s, np.abs(lhs-rhs).max(), (P*j).sum()))
