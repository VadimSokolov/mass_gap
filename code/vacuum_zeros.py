# Zeros in x = cos(theta) of the ground state of one plaquette, H = -lam*Laplacian + (1 - x)/lam, as an entire function of x.
# Zeros in a disc are counted by the argument principle and compared with the sign changes on the real axis below -1.
# High precision is needed: at 60 digits and lam = 0.1 the count on the circle is wrong.
import mpmath as mp
def coeffs(lam, N):
    H = mp.zeros(N, N)
    for i in range(N):
        H[i, i] = mp.mpf(lam)*((i+1)**2-1)/4+1/mp.mpf(lam)
        if i+1 < N: H[i, i+1] = H[i+1, i] = -1/(2*mp.mpf(lam))
    E, V = mp.eigsy(H); i0 = min(range(N), key=lambda i: E[i]); c = [V[k, i0] for k in range(N)]; return [-x for x in c] if c[0] < 0 else c
def psi(c, x):
    u0, u1 = mp.mpf(1), 2*x; s = c[0]*u0+c[1]*u1
    for k in range(2, len(c)): u0, u1 = u1, 2*x*u1-u0; s += c[k]*u1
    return s
def census(lam, R, N, M=6000, G=30000):
    c = coeffs(lam, N); tot = 0; prev = psi(c, mp.mpc(R, 0))
    for k in range(1, M+1): cur = psi(c, R*mp.expjpi(2*mp.mpf(k)/M)); tot += mp.arg(cur/prev); prev = cur
    xs = [-1-(R-1)*mp.mpf(k)/G for k in range(G+1)]; v = [psi(c, x) for x in xs]; z = [float(xs[i]) for i in range(G) if v[i]*v[i+1] < 0]; return round(float(tot/(2*mp.pi))), z
if __name__ == "__main__":
    mp.mp.dps = 160
    for lam, R, N in ((2.0, 12.0, 70), (0.5, 6.0, 70), (0.1, 2.5, 110)):
        w, z = census(lam, R, N); print("lambda = %g: zeros in |x| < %g: %d; real zeros below -1 in the disc: %d; nearest: %s" % (lam, R, w, len(z), [round(t, 4) for t in z[:4]]))
