# Exact rational certificate (standard library only).
# Karp-Sitnik, J. Approx. Theory 161 (2009), arXiv:math/0703084, Conjecture 2:
#   for 0 < sigma <= min a_i, sum (b_i - a_i) > 0, x > 0:
#   (q+1)F(q)(sigma, a; b; -x) < (1 + x K^(1/sigma))^(-sigma),  K = prod Gamma(a_i)Gamma(b_i - sigma)/(Gamma(b_i)Gamma(a_i - sigma)).
# Witness: q = 2, sigma = 1, a = (2, 2), b = (3/2, 3), x = 1/2.  At sigma = 1, Gamma(a)/Gamma(a-1) = a-1 and
# Gamma(b-1)/Gamma(b) = 1/(b-1), so K = prod (a_i-1)/(b_i-1) = (1/(1/2))*(1/2) = 1 and the bound is 2/3.
# 3F2(1,2,2; 3/2,3; -x) = sum_n t_n, t_n = (2)_n^2/((3/2)_n (3)_n) (-x)^n.  The ratio
# |t_{n+1}/t_n| = (n+2)^2 x/((n+3/2)(n+3)) <= 1/2 for every n >= 0 at x = 1/2, since (n+2)^2 <= (n+3/2)(n+3);
# so the terms alternate in sign and decrease geometrically in magnitude, and by the alternating-series bracket
# the sum lies between consecutive partial sums.  The loop only rechecks the first N ratios.
from fractions import Fraction as F
x = F(1, 2); N = 40
K = (F(2) - 1) / (F(3, 2) - 1) * (F(2) - 1) / (F(3) - 1)
assert K == 1
bound = 1 / (1 + x * K)
t = []; c = F(1)
for n in range(N + 1):
    t.append(c * (-x) ** n)
    c = c * (2 + n) ** 2 / ((F(3, 2) + n) * (3 + n))
for n in range(N):
    assert abs(t[n + 1]) <= abs(t[n]) / 2
lo = sum(t[:N + 1]) if t[N] < 0 else sum(t[:N])
hi = sum(t[:N]) if t[N] < 0 else sum(t[:N + 1])
assert lo > bound
print("3F2 in [%.12f, %.12f] > bound = %s (exact)" % (float(lo), float(hi), bound))
