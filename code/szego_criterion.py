# The Schur and Szego criterion for real zeros under convolution on SU(2).
# Theorem: for polynomials f, g in L, if f * f_D is in L' (D = deg g), then f * g is in L.
# The hypothesis is one statement about f alone: an antiderivative of f * f_D has only real zeros, all at most -1.
# Tested here for random and adversarial f; compared with the older hypothesis "f in L'"; and f * g checked directly.
import mpmath as mp, random
mp.mp.dps = 120
def poly_from_gam(gam):
    p=[mp.mpf(1)]
    for g in gam:
        q=[mp.mpf(0)]*(len(p)+1)
        for i,a in enumerate(p): q[i]+=a; q[i+1]+=a/mp.mpf(g)
        p=q
    return p  # ascending in x
def cheb_u(n):  # list of U_0..U_{n-1} ascending coeffs
    U=[[mp.mpf(1)],[mp.mpf(0),mp.mpf(2)]]
    for k in range(2,n):
        a=[mp.mpf(0)]+[2*t for t in U[-1]]; b=U[-2]+[mp.mpf(0)]*(len(a)-len(U[-2])); U.append([x-y for x,y in zip(a,b)])
    return U[:n]
def to_char(p):  # p ascending in x, degree d -> coefficients b_1..b_{d+1} with p = sum b_n U_{n-1}
    d=len(p)-1; U=cheb_u(d+1); r=p[:]; b=[mp.mpf(0)]*(d+1)
    for n in range(d,-1,-1):
        b[n]=r[n]/U[n][n]
        for i,t in enumerate(U[n]): r[i]-=b[n]*t
    return b
def from_char(b):
    U=cheb_u(len(b)); p=[mp.mpf(0)]*len(b)
    for n,c in enumerate(b):
        for i,t in enumerate(U[n]): p[i]+=c*t
    return p
def conv(p,q):
    a=to_char(p); b=to_char(q); m=min(len(a),len(b)); return from_char([a[n]*b[n]/(n+1) for n in range(m)])
def fpow(D): 
    p=[mp.mpf(1)]
    for _ in range(D): p=[mp.mpf(0)]+p; p=[p[i]+(p[i+1] if i+1<len(p) else 0) for i in range(len(p))]
    return [mp.binomial(D,k) for k in range(D+1)]
def roots(p):
    p=p[:]
    while len(p)>1 and abs(p[-1])<mp.mpf(10)**-60*max(abs(c) for c in p): p=p[:-1]
    return mp.polyroots(p[::-1],maxsteps=800,extraprec=600) if len(p)>1 else []
def pv(p,x): return sum(c*x**i for i,c in enumerate(p))
def in_L(p):
    r=roots(p); return all(abs(mp.im(z))<mp.mpf(10)**-25*(1+abs(z)) and mp.re(z)<=-1+mp.mpf(10)**-25 for z in r)
def in_Lprime(p):
    r=roots(p)
    if not all(abs(mp.im(z))<mp.mpf(10)**-25*(1+abs(z)) and mp.re(z)<=-1+mp.mpf(10)**-25 for z in r): return False
    pts=sorted([mp.re(z) for z in r]); F=[mp.mpf(0)]+[c/(i+1) for i,c in enumerate(p)]
    if p[-1]<0: F=[-c for c in F]
    vals=[pv(F,x) for x in pts]; D=len(pts)
    mins=[vals[D-1-i] for i in range(0,D,2)]; maxs=[vals[D-1-i] for i in range(1,D,2)]
    return max(mins)<=min(maxs+[pv(F,mp.mpf(-1))])+mp.mpf(10)**-40

if __name__ == "__main__":
    random.seed(3)
    print("%4s %4s %8s %12s %14s" % ("deg", "D", "laws", "f in L'", "f*f_D in L'"))
    for d, D in ((3, 3), (4, 4), (6, 6), (9, 9), (12, 12), (3, 6), (4, 12), (6, 9)):
        n = 120; a = b = 0
        for t in range(n):
            mode = t % 3
            if mode == 0: gam = [1+10**random.uniform(-2, 2) for _ in range(d)]
            elif mode == 1: gam = [1+10**random.uniform(-4, -1) for _ in range(d//2)]+[1+10**random.uniform(1, 3) for _ in range(d-d//2)]
            else: c = 1+10**random.uniform(-2, 2); e = 1+10**random.uniform(-2, 2); gam = [c]*(d//2)+[e]*(d-d//2)
            f = poly_from_gam(gam); a += in_Lprime(f); b += in_Lprime(conv(f, fpow(D)))
        print("%4d %4d %8d %12.2f %14.2f" % (d, D, n, a/n, b/n))
    random.seed(5); print("smoothing by a lower degree: g of degree D, f_d with d < D")
    for d, D in ((1, 7), (2, 6), (3, 9), (3, 12), (5, 8), (2, 12)):
        n = 80; b = 0
        for t in range(n):
            gam = [1+10**random.uniform(-2, 2) for _ in range(D)] if t % 2 else [1+10**random.uniform(-4, -1) for _ in range(D//2)]+[1+10**random.uniform(1, 3) for _ in range(D-D//2)]
            b += in_Lprime(conv(poly_from_gam(gam), fpow(d)))
        print("   D = %2d, d = %d: g*f_d in L' for %d of %d" % (D, d, b, n))
    random.seed(8); ok = tot = 0
    for df, D in ((2, 3), (3, 5), (2, 6), (3, 9), (5, 8), (1, 7), (4, 5), (7, 7)):
        for t in range(40):
            f = poly_from_gam([1+10**random.uniform(-2, 2) for _ in range(df)]); g = poly_from_gam([1+10**random.uniform(-2, 2) for _ in range(D)]); tot += 1; ok += in_L(conv(f, g))
    print("direct check of f * g in L: %d of %d pairs" % (ok, tot))
    # how large can D be before the hypothesis fails for an f outside L'?
    f = poly_from_gam([1.0, 10.0, 10.5]); print("f = (1+x)(1+x/10)(1+x/10.5) in L':", in_Lprime(f), "; f*f_D in L' for D =", [D for D in (3, 6, 12, 24, 48, 96) if in_Lprime(conv(f, fpow(D)))])
    print("bond moved laws rho^k, all zeros of even multiplicity: (zeros of rho, k, in L', f*f_d in L', f*f in L)")
    for gam, k in (([1.2, 2, 5], 4), ([1.05, 9], 2), ([1.01, 1.5, 30], 4), ([1.3, 1.31], 4), ([1.001, 1000], 4), ([2, 3, 4], 2)):
        f = poly_from_gam([c for c in gam for _ in range(k)]); d = len(gam)*k; print("  ", gam, k, in_Lprime(f), in_Lprime(conv(f, fpow(d))), in_L(conv(f, f)))
