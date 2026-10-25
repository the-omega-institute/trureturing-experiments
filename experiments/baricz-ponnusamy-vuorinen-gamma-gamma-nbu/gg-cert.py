# Rigorous rational certificate (standard library only), outward-rounded dyadic arithmetic.
# Baricz-Ponnusamy-Vuorinen, Expo. Math. 29 (2011), arXiv:1009.4814, Question 8: is the gamma-gamma distribution NBU,
# i.e. S(u1 + u2) <= S(u1) S(u2) for all u1, u2 > 0, with S the survival function?
# Witness a = b = 1, alpha = 1: density f(u) = 2 K_0(2 sqrt u), survival S(u) = int_0^oo exp(-s - u/s) ds.
# NBU at u1 = u2 = 1/2 needs S(1) <= S(1/2)^2; this script certifies S(1) > S(1/2)^2.
# g_u(s) = exp(-s - u/s) increases on (0, sqrt u] and decreases on [sqrt u, oo).  With rational split points
# pl <= sqrt u <= ph: on [s0, pl] (increasing) and [ph, T] (decreasing) Riemann sums bracket the integral;
# on [pl, ph] and on (0, s0] the integrand lies in [0, 1]; the tail beyond T is in [0, exp(-T)].
# exp(-y), y >= 0 rational, is enclosed as exp(-1)^m * exp(-r) with y = m + r, 0 <= r < 1, using a rational
# enclosure of exp(-1) and the alternating Taylor bracket for exp(-r); every product is rounded outward.
from fractions import Fraction as F
D = 1 << 96
def dn(v): return F(v.numerator * D // v.denominator, D)
def up(v): return F(-((-v.numerator * D) // v.denominator), D)
def mul(a, b):
    c = [a[0]*b[0], a[0]*b[1], a[1]*b[0], a[1]*b[1]]
    return (dn(min(c)), up(max(c)))
E1 = (F(36787944117, 10**11), F(36787944118, 10**11))   # exp(-1) = 0.367879441171...
def _alt_sum(n):
    # partial sum of sum_{j<=n} (-1)^j / j!; for n odd it is <= exp(-1), for n even it is >= exp(-1)
    s = F(0); t = F(1)
    for j in range(n + 1):
        s += t; t = -t / (j + 1)
    return s
assert E1[0] <= _alt_sum(23) and _alt_sum(22) <= E1[1], "exp(-1) enclosure not certified"
def exp_neg_r(r, n=24):
    # alternating series for exp(-r), 0 <= r < 1: partial sums alternate around the value
    s = F(0); t = F(1); sums = []
    for j in range(n):
        s += t; sums.append(s); t = t * (-r) / (j + 1)
    return (dn(min(sums[-1], sums[-2])), up(max(sums[-1], sums[-2])))
def exp_neg(y):
    m = y.numerator // y.denominator; r = y - m
    out = exp_neg_r(r)
    for _ in range(m): out = mul(out, E1)
    return out
def S_bounds(u, pl, ph, N=60, T=40):
    s0 = pl / N
    lo = hi = F(0)
    inc = [s0 + (pl - s0) * F(i, N) for i in range(N + 1)]
    dec = [ph + (T - ph) * F(i, 8 * N) for i in range(8 * N + 1)]
    g = lambda s: exp_neg(s + u / s)
    gi = [g(s) for s in inc]; gd = [g(s) for s in dec]
    for i in range(N):
        h = inc[i + 1] - inc[i]; lo += h * gi[i][0]; hi += h * gi[i + 1][1]
    for i in range(8 * N):
        h = dec[i + 1] - dec[i]; lo += h * gd[i + 1][0]; hi += h * gd[i][1]
    hi += s0 + (ph - pl)          # (0, s0] and [pl, ph], integrand <= 1
    hi += exp_neg(F(T))[1]        # tail
    return dn(lo), up(hi)
assert F(7071, 10**4)**2 <= F(1, 2) <= F(7072, 10**4)**2
lo1, hi1 = S_bounds(F(1), F(1), F(1))
lo5, hi5 = S_bounds(F(1, 2), F(7071, 10**4), F(7072, 10**4))
print('S(1)   in [%.6f, %.6f]' % (float(lo1), float(hi1)))
print('S(1/2) in [%.6f, %.6f]' % (float(lo5), float(hi5)))
assert lo1 > hi5 ** 2
print('certified: S(1) > S(1/2)^2, so NBU fails at u1 = u2 = 1/2 (exact rational comparison)')
