# Simulations for GBC matching. Mode "train": 4^4 lattices at couplings drawn from the prior, summaries (W11, W12, W22).
# Mode "block": 8^4 lattices at a given coupling, blocked to 4^4 by several rules, same summaries. Both average 16 configurations.
# usage: gbc_sim.py train SEEDIDX NSIM [REP]  |  gbc_sim.py block SEEDIDX BETA [REP]; REP (default 0 = original run) replicates; files go to $GBC_DIR (default /tmp)
import sys, os, numpy as np, su2_core as C
DIR = os.environ.get("GBC_DIR", "/tmp")
from mcrg_pl import block
def summ(U):
    pl = [(m, n) for m in range(4) for n in range(m+1, 4)]
    return [np.mean([C.loop(U, m, n, 1, 1)[..., 0].mean() for m, n in pl]),
            np.mean([C.loop(U, m, n, 1, 2)[..., 0].mean() for m in range(4) for n in range(4) if m != n]),
            np.mean([C.loop(U, m, n, 2, 2)[..., 0].mean() for m, n in pl])]
mode, seed = sys.argv[1], int(sys.argv[2]); SEED = int(sys.argv[4]) if len(sys.argv) > 4 else 0; C.rng = np.random.default_rng(1000+seed+100*SEED); out = []
if mode == "train":
    for b in C.rng.uniform(0.3, 9.0, int(sys.argv[3])):
        U = C.cold(4, 4)
        for i in range(40): U = C.sweep(U, b)
        acc = []
        for i in range(16):
            for j in range(3): U = C.sweep(U, b)
            acc.append(summ(U))
        out.append([b]+list(np.mean(acc, 0)))
    np.save(DIR+"/gbc_train_%d.npy" % seed, np.array(out))
else:
    b = float(sys.argv[3]); smears = [0.0, 0.1, 0.2, 0.3, 0.5]; U = C.cold(8, 4); acc = {s: [] for s in smears}
    for i in range(100): U = C.sweep(U, b)
    for i in range(16):
        for j in range(5): U = C.sweep(U, b)
        for s in smears: acc[s].append(summ(block(U, s)))
    np.save(DIR+"/gbc_block_%s.npy" % sys.argv[3], np.array([[b, s]+list(np.mean(acc[s], 0)) for s in smears]))
