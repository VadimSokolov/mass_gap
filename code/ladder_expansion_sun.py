# Weak coupling expansion of the projected ladder for SU(N), N = 2, 3, 4, in exact arithmetic.
# Small time heat kernel on the maximal torus: K_s ~ prod_{i<j} j(a_ij) exp(-|phi|^2/(2s)), j(x) = (x/2)/sin(x/2), a_ij = phi_i - phi_j,
# Haar ~ prod sin^2(a_ij/2). Time s is normalised by u_F = exp(-C_F s), C_F = (N^2-1)/(2N), and beta = N/s.
# After bond moving with power k the torus law is  prod a_ij^2 * prod j(a_ij)^(k-2) * exp(-|phi|^2/(2 tau)), tau = s/k.
import sympy as sp
def expansion(N):
    k, D, tau, s = sp.symbols("k D tau s", positive=True); m = k-2; z = sp.symbols("z1:%d" % N, real=True)
    basis = []
    for r in range(1, N):
        v = [sp.Integer(1)]*r+[-sp.Integer(r)]+[0]*(N-r-1); basis.append([x/sp.sqrt(r+r*r) for x in v])
    phi = [sp.expand(sum(z[i]*basis[i][j] for i in range(N-1))) for j in range(N)]
    roots = [phi[i]-phi[j] for i in range(N) for j in range(i+1, N)]
    L2 = sp.expand(sum(r**2 for r in roots)/24); L4 = sp.expand(sum(r**4 for r in roots)/2880); L6 = sp.expand(sum(r**6 for r in roots)/181440)
    P = lambda j: sp.expand(sum(p**j for p in phi))
    # graded pieces: index = power of tau
    w = [1, m*L2, m*L4+m**2*L2**2/2, m*L6+m**2*L2*L4+m**3*L2**3/6]
    c = [1, -P(2)/(2*N), P(4)/(24*N), -P(6)/(720*N)]
    vdm = sp.expand(sp.prod([r**2 for r in roots]))
    def gauss(poly):
        tot = 0
        for mon, co in sp.Poly(sp.expand(poly*vdm), *z).terms():
            if any(p % 2 for p in mon): continue
            tot += co*sp.prod([sp.factorial2(p-1) if p > 0 else 1 for p in mon])
        return sp.expand(tot)
    g0 = gauss(1)
    den = sum(gauss(w[i])/g0*tau**i for i in range(4))
    num = sum(gauss(sp.expand(sum(c[j]*w[i-j] for j in range(i+1))))/g0*tau**i for i in range(4))
    lnu = sp.series(sp.log(num)-sp.log(den), tau, 0, 4).removeO(); CF = sp.Rational(N*N-1, 2*N)
    return sp.expand((-D*lnu/CF).subs(tau, s/k)), (k, D, s)
for N in (2, 3, 4):
    e, (k, D, s) = expansion(N); co = sp.Poly(sp.expand(e*k/D), s).all_coeffs()[::-1]
    A, B = sp.factor(co[2]), sp.factor(co[3]); rho = sp.factor(sp.simplify((B-A**2)/A))
    print("SU(%d): s' = (D/k) s (1 + A s + B s^2),  A = %s,  B = %s,  two loop ratio (B - A^2)/A = %s" % (N, A, B, rho), flush=True)
