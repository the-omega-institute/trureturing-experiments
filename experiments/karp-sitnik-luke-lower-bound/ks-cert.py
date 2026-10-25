# Exact rational certificate (standard library only).
# Karp-Sitnik, J. Approx. Theory 161 (2009), arXiv:math/0703084, Conjecture 1:
# "Theorem 3 is true for all sigma > 0 and sum_i (b_i - a_i) > 0", Theorem 3 being
#   1/(1 + x prod a_i/b_i)^sigma < (q+1)F(q)(sigma, a; b; -x).
# Witness: q = 2, sigma = 1, a = (1, 1), b = (1/10, 2), x = 1/10; sum(b - a) = 1/10 > 0, prod a/b = 5.
# 3F2(1,1,1; 1/10, 2; -x) = sum_n t_n with t_n = n!/((1/10)_n (n+1)) (-x)^n.
# For n >= 1 the ratio |t_{n+1}/t_n| = (n+1)/(n+1/10) * (n+1)/(n+2) * x is < 1 and decreasing to x,
# so the series is alternating with eventually decreasing terms; partial sums bracket the value.
from fractions import Fraction as F
x = F(1, 10); N = 40
t = []; c = F(1)
for n in range(N + 1):
    t.append(c / (n + 1) * (-x) ** n)
    c = c * (1 + n) / (F(1, 10) + n)
# check the alternating/decreasing property on the tail start
for n in range(1, N):
    assert abs(t[n + 1]) < abs(t[n])
S_lo = sum(t[:N + 1]) if t[N] < 0 else sum(t[:N])   # ends on a negative term -> lower bound
S_hi = sum(t[:N]) if t[N] < 0 else sum(t[:N + 1])
bound = 1 / (1 + x * F(5))
assert S_hi < bound
print("3F2 in [%.12f, %.12f] < bound = %s (exact)" % (float(S_lo), float(S_hi), bound))
