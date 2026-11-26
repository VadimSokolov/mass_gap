# An inequality of type (M) for the ladder, tested over two blockings: SU(2) on 16^4, blocked to 8^4 and to 4^4 (scale factor 2, smear 0.5),
# count of the plaquette at each level against n exact rungs of the ladder from the Wilson law matched to the fine plaquette.
# Counts from characters, 2 sum_{n<=40} E chi_n - 2. Usage: beta [configs]
import sys, numpy as np
from scipy.special import ive
from scipy.optimize import brentq
import su2_core as C
from mcrg_pl import block
from count_blocked import char_count, plaq
TH = (np.arange(6000)+0.5)*np.pi/6000; HA = np.sin(TH)**2; NN = np.arange(1, 241); CH = np.sin(np.outer(NN, TH))/np.sin(TH)[None, :]
def ladder(kappa, rungs):
    rho = np.exp(kappa*(np.cos(TH)-1)); out = []
    for r in range(rungs+1):
        u = (CH*(rho*HA)[None, :]).sum(1)/((rho*HA).sum()*NN); out.append((2*(NN*u).sum()-2, u[1]))
        b = rho**4; ub = (CH*(b*HA)[None, :]).sum(1)/((b*HA).sum()*NN); rho = np.maximum((NN*ub**4)@CH, 1e-300)
    return out
if __name__ == "__main__":
    beta = float(sys.argv[1]); nc = int(sys.argv[2]) if len(sys.argv) > 2 else 30; C.rng = np.random.default_rng(int(100*beta)+3); U = C.cold(16, 4)
    for i in range(150): U = C.sweep(U, beta)
    A = [[], [], []]
    for i in range(nc):
        for q in range(5): U = C.sweep(U, beta)
        B1 = block(U, 0.5); B2 = block(B1, 0.5)
        for lev, V in enumerate((U, B1, B2)): A[lev].append(plaq(V))
    res = []
    for v in A:
        v = np.array(v); c = np.array([char_count(v[g::5].ravel()) for g in range(5)]); res.append((char_count(v.ravel()), c.std()/np.sqrt(5), v.mean()))
    kap = brentq(lambda k: ive(2, k)/ive(1, k)-res[0][2], 1e-6, 400); lad = ladder(kap, 2)
    print("beta %.2f  " % beta+"  ".join("level %d: lattice %.3f (%.3f) u2 %.4f | ladder %.3f u2 %.4f" % (n, res[n][0], res[n][1], res[n][2], lad[n][0], lad[n][1]) for n in range(3)))
