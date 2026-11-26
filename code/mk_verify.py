"""Heat kernel projected MK map for SU(2): R(beta), robustness, and grid verification of R(beta) < beta.
Coordinate: beta = 2/s with u_2 = exp(-3s/4) (heat kernel). Polson and Sokolov, Nov 2026."""
import numpy as np, sys
def setup(M, NMAX):
    x, w = np.polynomial.legendre.leggauss(M)
    th = (x+1)*np.pi/2; w = w*np.pi/2; haar = (2/np.pi)*np.sin(th)**2*w
    n = np.arange(1, NMAX+1); CHI = np.sin(np.outer(n, th))/np.sin(th)
    return th, haar, n, CHI
def make_R(M=6000, NMAX=1500, b=2):
    th, haar, n, CHI = setup(M, NMAX); k = b**2; b2 = b**2
    def R(beta):
        s = 2/beta; u = np.exp(-s*(n**2-1)/4)
        f = np.clip((n*u)@CHI, 1e-300, None); lf = k*np.log(f)
        g = np.exp(lf-lf.max()); c = (CHI@(g*haar))/n; c = c/c[0]
        logu2 = b2*np.log(c[1])
        return 2/(-(4/3)*logu2)            # back to heat kernel beta
    return R
if __name__ == "__main__":
    R = make_R()
    for bb in [0.05, 0.2, 0.5, 1, 2, 5, 10, 20, 50]:
        print("beta=%6.2f  R=%.6f  delta=%.6f" % (bb, R(bb), bb-R(bb)))
    # robustness to truncation
    R2 = make_R(4000, 800)
    print("truncation check max |dR|:", max(abs(R(b)-R2(b)) for b in [0.1, 1, 5, 20, 50]))
    # b = 3: predicted 8/27
    R3 = make_R(6000, 1500, b=3)
    for bb in [10, 20, 40]: print("b=3 beta=%g delta=%.5f (8/27=%.5f, AF 11ln3/3pi^2=%.5f)" % (bb, bb-R3(bb), 8/27, 11*np.log(3)/(3*np.pi**2)))
