# Absolute string tension from the Migdal Kadanoff ladder, Wilson start, SU(2), d = 3 and 4.
# sigma a^2 = lim_m b^{-2m} (-log u_2^{(m)}), with u_n the normalised character coefficients after m steps.
import numpy as np
from scipy.special import ive
M, NMAX = 4000, 400
x, w = np.polynomial.legendre.leggauss(M); th = (x+1)*np.pi/2; w = w*np.pi/2
haar = (2/np.pi)*np.sin(th)**2*w; n = np.arange(1, NMAX+1); CHI = np.sin(np.outer(n, th))/np.sin(th)
def sigma(beta, d=4, b=2, steps=40, stop=25.0):
    # stop: leave the loop once -log u_2 exceeds this. One more rung multiplies -log u_2 by about b^2, so the last c_2/c_1 is at least about exp(-b^2*stop) = 1e-11 for stop = 25, still far above the double precision noise (1e-16).
    # With the old value 60 a rung starting at -log u_2 = 49 gave c_2/c_1 ~ exp(-195), below double precision: the 1e-300 clip then gave 2760 and the result jumped by about 4.
    logu = np.log(np.clip(ive(n, beta)/ive(1, beta), 1e-300, None)); k = b**(d-2); scale = 1.0
    for m in range(steps):
        f = np.clip((n*np.exp(logu))@CHI, 1e-300, None); lf = k*np.log(f); g = np.exp(lf-lf.max())
        c = (CHI@(g*haar))/n; c = np.clip(c/c[0], 1e-300, None)
        logu = b*b*np.log(c); scale *= b*b
        if -logu[1] > stop: break
    assert -logu[1] < 150, "u_2 reached the double precision floor: lower stop"
    return -logu[1]/scale
if __name__ == "__main__":
    data4 = {2.2: 0.4690, 2.3: 0.3690, 2.4: 0.2660, 2.5: 0.1871, 2.635: 0.1207}
    print("d=4: beta, ladder a sqrt(sigma), lattice, ratio")
    for be, v in data4.items():
        s = np.sqrt(sigma(be)); print("  %.3f  %.4f  %.4f  %.2f" % (be, s, v, s/v))
    print("d=3: beta, ladder beta a sqrt(sigma) (lattice continuum value 1.341)")
    for be in [4, 6, 9, 12, 18, 24]:
        print("  %5.1f  %.4f" % (be, be*np.sqrt(sigma(be, d=3, steps=80))))
    np.save("/tmp/mk_sigma4.npy", np.array([[be, np.sqrt(sigma(be))] for be in np.linspace(0.5, 3.2, 28)]))
    np.save("/tmp/mk_sigma3.npy", np.array([[be, np.sqrt(sigma(be, d=3, steps=80))] for be in np.linspace(1, 24, 24)]))
