# Outward-rounded rational certificate (standard library only).
# Kalmykov-Karp, arXiv:1702.03788, Conjecture 1, coefficient clause.
# g_n = 1, q = r^10 with r = 99/100, mu = -7/10, alpha = beta = 7/10, coefficient of x^2:
#   delta_2 = 1 - P*C,
#   P = (q^{3/10};q)_inf (q^{7/10};q)_inf / ((q;q)_inf^2 (1-q)),
#   C = [-7/10]_q/([7/10]_q [17/10]_q) + 1/[7/10]_q + 1/[3/10]_q,  [x]_q = (1-q^x)/(1-q).
# Infinite products: N = 800 factors, rounded outward to multiples of 2^-200,
# with prod_{j>=N} (1 - a q^j) >= 1 - a q^N/(1-q) for the lower bound.
from fractions import Fraction as F
r = F(99, 100); q = r**10; D = 1 << 200
def dn(x): return F(x.numerator * D // x.denominator, D)
def up(x): return F(-((-x.numerator * D) // x.denominator), D)
def br(n): return (1 - r**n) / (1 - q)          # [n/10]_q
C = br(-7) / (br(7) * br(17)) + 1 / br(7) + 1 / br(3)
def bounds(a, N=800):
    lo = F(1); hi = F(1); xl = a; xh = a
    for _ in range(N):
        lo = dn(lo * (1 - xh)); hi = up(hi * (1 - xl)); xl = dn(xl * q); xh = up(xh * q)
    return dn(lo * (1 - xh / (1 - q))), hi
l3, _ = bounds(r**3); l7, _ = bounds(r**7); _, uq = bounds(q)
Plow = l3 * l7 / (uq * uq * (1 - q))
assert C > 0
assert Plow * C > F(1064, 1000)
print("C > 0 and P*C > 1.064 (exact); delta_2 <=", float(1 - Plow * C), "(approx)")
