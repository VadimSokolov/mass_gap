# A first empirical rung: block 8^4 configurations to 4^4 and fit a single plaquette coupling by pseudolikelihood.
# Also the smallest Griffiths test, W_12 against W_11^2. Usage: beta
import sys, numpy as np
from scipy.special import ive
from scipy.optimize import brentq
import su2_core as C
def mple(U):
    D = U.shape[-2]; ks, a0 = [], []
    for mu in range(D):
        A = C.staple(U, mu); k = np.linalg.norm(A, axis=-1); ks.append(k.ravel()); a0.append(C.qm(U[..., mu, :], A/k[..., None])[..., 0].ravel())
    ks, a0 = np.concatenate(ks), np.concatenate(a0)
    return brentq(lambda b: np.sum(ks*(a0-ive(2, b*ks)/ive(1, b*ks))), 0.01, 200)
def block(U, smear):
    """scale factor 2: blocked link = product of two links, plus smear times the sum of the 2 x 1 staples, projected to SU(2)"""
    D = U.shape[-2]; out = []
    for mu in range(D):
        V = C.line(U, mu, 2); S = 0
        for nu in range(D):
            if nu == mu: continue
            Un = U[..., nu, :]; up = C.qm(C.qm(Un, C.sh(C.line(U, mu, 2), nu, 1)), C.cj(C.sh(C.sh(Un, mu, 1), mu, 1)))
            Und = C.sh(Un, nu, -1); dn = C.qm(C.qm(C.cj(Und), C.sh(C.line(U, mu, 2), nu, -1)), C.sh(C.sh(Und, mu, 1), mu, 1))
            S = S + up + dn
        V = V + smear*S; V = V/np.linalg.norm(V, axis=-1)[..., None]
        out.append(V[(slice(None, None, 2),)*D])
    return np.stack(out, axis=-2)
if __name__ == "__main__":
    beta = float(sys.argv[1]); U = C.cold(8, 4)
    for i in range(100): U = C.sweep(U, beta)
    f, b0, b1, w11, w12 = [], [], [], [], []
    for i in range(40):
        for j in range(5): U = C.sweep(U, beta)
        f.append(mple(U)); b0.append(mple(block(U, 0.0))); b1.append(mple(block(U, 0.5)))
        w11.append(np.mean([C.loop(U, m, n, 1, 1)[..., 0].mean() for m in range(4) for n in range(m+1, 4)]))
        w12.append(np.mean([C.loop(U, m, n, 1, 2)[..., 0].mean() for m in range(4) for n in range(4) if m != n]))
    se = lambda x: np.std(x)/np.sqrt(len(x))
    print("beta %.2f  MPLE fine %.3f(%.3f)  blocked plain %.3f(%.3f)  blocked smeared %.3f(%.3f)  W12 %.5f  W11^2 %.5f  diff %.5f(%.5f)" % (beta, np.mean(f), se(f), np.mean(b0), se(b0), np.mean(b1), se(b1), np.mean(w12), np.mean(w11)**2, np.mean(w12)-np.mean(w11)**2, se(np.array(w12)-np.array(w11)**2)), flush=True)
