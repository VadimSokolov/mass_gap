# The count of the plaquette of a blocked lattice (8^4 to 4^4, scale factor 2, with and without smearing), against one rung of the exact ladder.
# No link can be integrated out of a smeared plaquette, so the count is taken from the characters: count + 2 = 2 sum_n E[chi_n]/1, truncated at n = 40
# (each term has finite variance; the plain average of 1/(1 - tr U/2) has not). Checked on the unsmeared case against count_lattice.py. Usage: beta [configs]
import sys, numpy as np
from scipy.special import ive, gammaln
from scipy.optimize import brentq
import su2_core as C
from mcrg_pl import block
from count_lattice import wilson_count, heat_count
from count_lyapunov import rung, logZ, j, pad
NN = 40
def char_count(a0):
    th = np.arccos(np.clip(a0, -1+1e-12, 1-1e-12)); return 2*sum((np.sin(n*th)/np.sin(th)).mean() for n in range(1, NN+1))-2
def plaq(U): return np.concatenate([C.loop(U, m, n, 1, 1)[..., 0].ravel() for m in range(4) for n in range(m+1, 4)])
if __name__ == "__main__":
    beta = float(sys.argv[1]); nc = int(sys.argv[2]) if len(sys.argv) > 2 else 40; C.rng = np.random.default_rng(int(100*beta)+7); U = C.cold(8, 4)
    for i in range(80): U = C.sweep(U, beta)
    A = {s: [] for s in ("fine", 0.0, 0.5, 2.0)}
    for i in range(nc):
        for q in range(4): U = C.sweep(U, beta)
        A["fine"].append(plaq(U))
        for s in (0.0, 0.5, 2.0): A[s].append(plaq(block(U, s)))
    res = {}
    for s, v in A.items():
        v = np.array(v); c = np.array([char_count(v[g::5].ravel()) for g in range(5)]); res[s] = (char_count(v.ravel()), c.std()/np.sqrt(5), v.mean())
    kap = brentq(lambda k: ive(2, k)/ive(1, k)-res["fine"][2], 1e-6, 400); l = j*np.log(kap)-gammaln(j+1)+logZ; m0, mb, m1, u2 = rung(pad(np.exp(l-l.max())))
    print("beta %.2f fine: count %.3f (%.3f) u2 %.4f | ladder from the matched Wilson law: count %.3f -> %.3f, u2 -> %.4f" % (beta, res["fine"][0], res["fine"][1], res["fine"][2], m0, m1, u2))
    for s in (0.0, 0.5, 2.0): print("   blocked, smear %.1f: count %.3f (%.3f)  u2 %.4f  Wilson matched %.3f  heat kernel matched %.3f" % (s, res[s][0], res[s][1], res[s][2], wilson_count(res[s][2]), heat_count(res[s][2])))
