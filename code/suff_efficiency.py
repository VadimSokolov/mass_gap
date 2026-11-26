# Sufficiency and efficiency in two dimensional SU(2): three estimators of beta from one configuration of a 16 x 16 torus.
#   MLE: solve u(beta) = mean plaquette (the sufficient statistic);  MPLE: Besag's pseudolikelihood on links;
#   lossy: solve P_beta(a0 > 1/2) = observed fraction (a one bit summary of each plaquette). Cramer Rao bound for comparison.
import numpy as np, sys
sys.path.insert(0, "code"); import su2_core as C
from scipy.special import ive
from scipy.optimize import brentq
from scipy.integrate import quad
u = lambda b: ive(2, b)/ive(1, b); info = lambda b: 1-3*u(b)/b-u(b)**2
def frac_half(b):
    f = lambda a: np.sqrt(1-a*a)*np.exp(b*(a-1)); return quad(f, 0.5, 1)[0]/quad(f, -1, 1)[0]
def mple(U):
    ks, a0 = [], []
    for mu in range(2):
        A = C.staple(U, mu); k = np.linalg.norm(A, axis=-1); ks.append(k.ravel()); a0.append(C.qm(U[..., mu, :], A/k[..., None])[..., 0].ravel())
    ks, a0 = np.concatenate(ks), np.concatenate(a0); return brentq(lambda b: np.sum(ks*(a0-u(b*ks))), 0.05, 60)
L = 16; NREP = 6
def chain(args):
    rep, beta = args
    C.rng = np.random.default_rng(int(beta*7)+1000*rep); U = C.cold(L, 2); est = []     # rep 0 is the configuration of the first version
    for i in range(50): U = C.sweep(U, beta)
    for i in range(300):
        for j in range(6): U = C.sweep(U, beta)
        p = C.loop(U, 0, 1, 1, 1)[..., 0]
        est.append((brentq(lambda b: u(b)-p.mean(), 0.01, 80), mple(U), brentq(lambda b: frac_half(b)-np.clip((p > 0.5).mean(), 0.02, 0.98), 0.01, 80)))
    est = np.array(est); cr = 1/np.sqrt(L*L*info(beta)); return rep, beta, est.mean(0), est.std(0, ddof=1), cr
if __name__ == "__main__":
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor() as ex: res = list(ex.map(chain, [(r, b) for r in range(NREP) for b in (1.0, 2.0, 4.0)]))
    for beta in (1.0, 2.0, 4.0):
        rr = [x for x in res if x[1] == beta]; cr = rr[0][4]; r0 = [x for x in rr if x[0] == 0][0]
        print("beta %.0f  replicate 0 (300 configurations): mean (MLE, MPLE, lossy) %s  sd %s  Cramer Rao %.4f  efficiency %s" % (beta, np.round(r0[2], 3), np.round(r0[3], 4), cr, np.round((cr/r0[3])**2, 2)))
        sd = np.array([x[3] for x in rr]); eff = (cr/sd)**2; se = eff.std(0, ddof=1)/np.sqrt(len(rr)); sdm = sd.mean(0); sdse = sd.std(0, ddof=1)/np.sqrt(len(rr))
        print("beta %.0f  over %d replicates (300 configurations each): sd %s (s.e. %s)  efficiency %s (s.e. %s)" % (beta, len(rr), np.round(sdm, 4), np.round(sdse, 4), np.round(eff.mean(0), 3), np.round(se, 3)), flush=True)
