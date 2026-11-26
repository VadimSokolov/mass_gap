# Exact zero free sector for two dimensional SU(2) with the heat kernel action on a torus of A plaquettes:
#   Z_A(t) = sum_{n>=1} exp(-(n^2-1) A t/4) = g(q),  q = exp(-A t/4),  g(q) = 1 + q^3 + q^8 + q^15 + ... = (theta_3(q) - 1)/(2q).
# Zeros of Z_A in t are t = 4 s_j/A with s_j = -log q_j, q_j the zeros of g in the unit disc: the same directions for every A.
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from mpmath import mp, mpc, findroot
mp.dps = 30
# g(q) = sum_{n>=1} q^(n^2-1), evaluated as a power series (valid for every |q|<1; 150 terms give |q|^(n^2) < 1e-29 for |q| <= 0.997).
NN = np.arange(1, 151)
def gs(q):
    q = np.asarray(q, complex)[..., None]; return (q**(NN*NN-1)).sum(-1)
def dgs(q):
    q = np.asarray(q, complex)[..., None]; return ((NN*NN-1)*q**(NN*NN-2)).sum(-1)
RMAX = 0.997
def winding(r, N=2**21):
    """Number of zeros of g in |q| < r (argument principle on the circle |q| = r). g(r e^{i th}) on N equispaced points by one FFT of its power series."""
    c = np.zeros(N); c[NN*NN-1] = r**(NN*NN-1.0); gv = np.fft.ifft(c)*N
    ph = np.unwrap(np.angle(np.append(gv, gv[0]))); assert np.abs(np.diff(ph)).max() < 2.0, "phase step not resolved on |q|=%.4f" % r
    return int(round((ph[-1]-ph[0])/(2*np.pi)))
# Step 1: argument principle counts on circles (an independent count that does not use any root finder).
radii = [0.80, 0.83, 0.85, 0.86, 0.87, 0.88, 0.90, 0.92, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995, RMAX]
counts = {r: winding(r) for r in radii}
print("argument principle, number of zeros of g in |q| < r (both half planes):", counts)
# Step 2: dense polar grid of starting points in the upper half disc, damped Newton on the series; failures are counted, not hidden.
R0, A0 = np.meshgrid(np.concatenate([np.linspace(0.5, 0.98, 60), np.linspace(0.981, RMAX, 80)]), np.linspace(0.003, np.pi, 800)); Q = (R0*np.exp(1j*A0)).ravel(); conv = np.zeros(Q.shape, bool); act = np.arange(Q.size)
for it in range(80):                    # damped Newton on the active starts only; converged or escaped starts are removed
    if act.size == 0: break
    q = Q[act]; d = dgs(q); step = gs(q)/np.where(d != 0, d, 1); big = np.abs(step) > 0.05; step[big] *= 0.05/np.abs(step[big]); qn = q - step
    Q[act] = np.where(np.abs(qn) < 0.9999, qn, q); done = np.abs(step) < 1e-14; esc = np.abs(qn) >= 0.9999
    conv[act[done & ~esc]] = True; act = act[~done & ~esc]
conv &= np.abs(gs(Q)) < 1e-10
print("Newton starts: %d, converged to a zero %d, left the disc or did not converge (dropped, counted here) %d" % (Q.size, conv.sum(), (~conv).sum()))
cand = []
for q in sorted(Q[conv], key=lambda z: (abs(z), z.imag)):
    if q.imag > 1e-9 and abs(q) <= RMAX+1e-9 and all(abs(q-w) > 1e-7 for w in cand): cand.append(q)
zs = []
for q in cand:
    try: zs.append(complex(findroot(lambda z: sum(z**(int(n)**2-1) for n in range(1, 151)), mpc(q), tol=1e-22)))
    except Exception as e: print("polish failed at", q, e); zs.append(q)
# Step 3: consistency, each disc must contain exactly (argument principle count)/2 upper half plane zeros (g has no real zeros: checked below).
xr = np.linspace(-0.999, 0.999, 200001); assert gs(xr).real.min() > 0.1 and np.abs(gs(xr).imag).max() < 1e-12, "real zero or non real value on the real axis"
for r in radii:
    n_found = 2*sum(1 for z in zs if abs(z) < r)
    print("  |q| < %.3f: argument principle %3d, Newton list %3d %s" % (r, counts[r], n_found, "" if n_found == counts[r] else "  MISMATCH"))
    assert n_found == counts[r], "root search incomplete below |q| = %.3f" % r
zs = np.array(sorted(zs, key=abs)); s = -np.log(zs); ang = np.degrees(np.abs(np.angle(s)))
for z, ss, an in zip(zs, s, ang): print("q = %.5f%+.5fi  |q| = %.5f   s = %.4f%+.4fi   angle of s from the real axis %.2f deg" % (z.real, z.imag, abs(z), ss.real, ss.imag, an))
print("smallest angle: %.2f degrees, from the zero with |q| = %.5f" % (ang.min(), abs(zs[ang.argmin()])))
# What the finite search must cover.  In the closed sector |arg s| <= phi, for |s| <= s0 the modular identity gives
#   |theta_3(e^-s)| >= sqrt(pi/|s|) (1 - 2 sum_n exp(-n^2 pi^2 cos(phi)/|s|)) > 1,  so g has no zero there;
# for |s| >= s0 the sector lies in |q| = exp(-Re s) <= exp(-s0 cos(phi)).  So the zeros that matter all lie in |q| <= RNEED below.
phi0 = np.radians(ang.min())
for s0 in np.arange(0.3, 1.5, 0.05):
    lb = np.sqrt(np.pi/s0)*(1-2*sum(np.exp(-k*k*np.pi**2*np.cos(phi0)/s0) for k in range(1, 60)))
    if lb > 1: best = (s0, lb)
s0, lb = best; RNEED = np.exp(-s0*np.cos(phi0))
print("modular bound: |theta_3| >= %.3f > 1 for |s| <= %.2f in the sector of half angle %.2f deg; so only zeros with |q| <= %.4f matter, and %d zeros of g (both half planes) lie there" % (lb, s0, ang.min(), RNEED, 2*sum(1 for z in zs if abs(z) <= RNEED)))
from scipy.optimize import brentq
print("Rouche radius: |g-1| <= r^3/(1-r^5) < 1 for r < %.4f (so 0.83 is provable by hand, 0.86 is not); true largest zero free radius = smallest |q_j| = %.5f" % (brentq(lambda r: r**3/(1-r**5)-1, 0.5, 0.99), abs(zs[0])))
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.3)); t = np.linspace(0, 2*np.pi, 300)
ax[0].plot(np.cos(t), np.sin(t), "k", lw=0.6); ax[0].plot(zs.real, zs.imag, "o", color="#c0392b", ms=4); ax[0].plot(zs.real, -zs.imag, "o", color="#c0392b", ms=4)
ax[0].set_aspect("equal"); ax[0].set_xlabel("Re $q$"); ax[0].set_ylabel("Im $q$"); ax[0].set_title("zeros of $1+q^3+q^8+\\cdots$", fontsize=9)
phi = np.radians(ang.min()); rr = np.linspace(0, 6, 2)
ax[1].fill_between(rr, -rr*np.tan(phi), rr*np.tan(phi), color="#aed6f1", label="zero free sector")
for A, c in [(4, "#1f4e79"), (9, "#16a085"), (25, "#8e44ad")]:
    b = 2/(4*s/A); bb = np.concatenate([b, np.conj(b)]); ax[1].plot(bb.real, bb.imag, "o", ms=3, color=c, label=f"zeros, $A={A}$")
ax[1].set_xlim(0, 6); ax[1].set_ylim(-6, 6); ax[1].set_xlabel(r"Re $\beta$"); ax[1].set_ylabel(r"Im $\beta$"); ax[1].legend(frameon=False, fontsize=7)
fig.tight_layout(); fig.savefig("figures/torus_sector.pdf")
