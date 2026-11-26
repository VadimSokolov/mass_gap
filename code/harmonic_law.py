# The harmonic law for the mean count of f_{a_1} * ... * f_{a_p}, in exact rational arithmetic.
# Exact: 2 sum_n n prod_i u_n(J_i) - 2, with J_i = a_i + 1.  Law: 1/A - 7/6 + B/(2A^2) - C/(3A^3), A, B, C = sum of J_i^(-1), J_i^(-2), J_i^(-3).
from fractions import Fraction as F
from math import factorial
import random
def un(n, J): return F(factorial(J-1)*factorial(J+1), factorial(J-n)*factorial(J+n)) if n <= J else F(0)
def exact(Js):
    tot = F(0)
    for n in range(1, min(Js)+1):
        t = F(n)
        for J in Js: t *= un(n, J)
        tot += t
    return 2*tot-2
def law(Js):
    A = sum(F(1, J) for J in Js); B = sum(F(1, J*J) for J in Js); C = sum(F(1, J**3) for J in Js); return 1/A-F(7, 6)+B/(2*A*A)-C/(3*A**3)
if __name__ == "__main__":
    random.seed(1); worst = 0
    for scale in (20, 80, 320):
        for p in (2, 3, 4, 8):
            w = 0
            for trial in range(12):
                Js = [random.randint(scale, 3*scale) for _ in range(p)]; w = max(w, abs(float(exact(Js)-law(Js)))*min(Js))
            print("scale %4d  p = %d  largest |exact - law| times the smallest J: %.4f" % (scale, p, w)); worst = max(worst, w)
    print("point masses, k = p = 4: drift", [round(float(exact([4*m+1]*4)-m), 5) for m in (20, 80, 320)], "limit -13/16 =", -13/16)
    print("largest scaled error", round(worst, 4))
