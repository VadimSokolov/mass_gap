# W12 - W11^2 at strong coupling in four dimensional SU(2), against the leading strong coupling term 2(d-2) u^8, u = W11.
import sys, numpy as np, su2_core as C
beta = float(sys.argv[1]); C.rng = np.random.default_rng(int(10*beta)); U = C.cold(8, 4); w11, w12 = [], []
for i in range(60): U = C.sweep(U, beta)
for i in range(400):
    for j in range(2): U = C.sweep(U, beta)
    w11.append(np.mean([C.loop(U, m, n, 1, 1)[..., 0].mean() for m in range(4) for n in range(m+1, 4)]))
    w12.append(np.mean([C.loop(U, m, n, 1, 2)[..., 0].mean() for m in range(4) for n in range(4) if m != n]))
w11, w12 = np.array(w11), np.array(w12); u = w11.mean(); bl = (w12-2*u*w11+u*u).reshape(20, 20).mean(1)   # linearised difference, blocked
print("beta %.1f  u %.4f  W12-W11^2 %.6f +- %.6f   4u^8 %.6f" % (beta, u, w12.mean()-u*u, bl.std()/np.sqrt(20), 4*u**8), flush=True)
