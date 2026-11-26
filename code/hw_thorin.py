# The Hartman Watson law is a generalized gamma convolution: its Thorin measure is -theta'(y) dy / pi at t = y^2/2,
# with theta(y) = arg I_{iy}(r) and theta'(y) = -(2/pi) int_0^{pi/2} K_0(2 r cos a) cosh(2 y a) da / |I_{iy}(r)|^2  (Neumann's product formula).
import numpy as np, mpmath as mp
from scipy.integrate import quad
from scipy.special import k0
mp.mp.dps = 30
def Iiy(y, r): return complex(mp.besseli(mp.mpc(0, y), r))
def dtheta(y, r):
    f = lambda a: k0(2*r*np.cos(a))*np.cosh(2*y*a); return -(2/np.pi)*quad(f, 0, np.pi/2, limit=200)[0]/abs(Iiy(y, r))**2
def phase(r, ys):
    out, prev = [], 0.0
    for y in ys:
        a = float(np.angle(Iiy(y, r))) if y > 0 else 0.0
        a += 2*np.pi*np.round((prev-a)/(2*np.pi)); out.append(a); prev = a
    return np.array(out)
def log_ratio(r, shift, s, Y=60.0, n=6000):
    # -int log(1 + s/(t + shift)) U(dt), t = y^2/2, with the tail beyond Y from theta' ~ -log(2y/r)
    y = (np.arange(n)+0.5)*Y/n; d = np.array([dtheta(v, r) for v in y]); return (np.log1p(2*s/(y*y+2*shift))*d).sum()*Y/n/np.pi-2*s*(np.log(2*Y/r)+1)/Y/np.pi
if __name__ == "__main__":
    for r in (0.5, 2.0, 8.0):
        num = [(float(np.angle(Iiy(y+1e-4, r)/Iiy(y-1e-4, r)))/2e-4, dtheta(y, r)) for y in (0.5, 2.0, 7.0)]
        print("r = %g  theta'(y), numerical and by the formula:" % r, ["%.6f %.6f" % p for p in num], " theta'(0) = -K0/I0 = %.6f" % float(-mp.besselk(0, r)/mp.besseli(0, r)))
    for r, s in ((1.0, 0.5), (1.0, 4.0), (4.0, 2.0)):
        print("r = %g s = %g  log(I_sqrt(2s)/I_0) = %.5f   from the Thorin measure %.5f" % (r, s, float(mp.log(mp.besseli(mp.sqrt(2*s), r)/mp.besseli(0, r))), log_ratio(r, 0.0, s)))
    for kap, nn in ((2.0, 2), (2.0, 4), (6.0, 3)):
        print("Wilson kappa = %g: log(I_%d/I_1) = %.5f   from the shifted Thorin measure %.5f" % (kap, nn, float(mp.log(mp.besseli(nn, kap)/mp.besseli(1, kap))), log_ratio(kap, 0.5, (nn*nn-1)/2)))
    r = 1.0; ys = np.linspace(0, 14, 5601); th = phase(r, ys); kz = [float(mp.findroot(lambda v: mp.besselk(mp.mpc(0, v), r).real, g)) for g in (2.96, 4.53, 5.88, 7.1, 8.3, 9.4, 10.5, 11.5, 12.5, 13.5)]
    kz = sorted(set(round(z, 4) for z in kz if 0 < z < 14)); print("zeros in y of K_iy(1):", kz); print("-theta/pi at those zeros:  ", [round(float(-np.interp(z, ys, th)/np.pi), 4) for z in kz])
