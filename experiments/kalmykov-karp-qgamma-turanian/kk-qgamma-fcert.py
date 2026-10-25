# Outward-rounded rational interval certificate (standard library only).
# Kalmykov-Karp, arXiv:1702.03788, Theorem 3 / Conjecture 1, function clause:
#   Delta_g(alpha,beta;x) = g(mu+alpha;x) g(mu+beta;x) - g(mu;x) g(mu+alpha+beta;x) > 0 for 0 < x < R_g.
# Witness: g_n = 1, q = r^10 with r = 93/100, mu = -7/10, alpha = 3/10, beta = 7/10, x = X (default 93/100).
# Normalisation: h(s;x) = (q;q)_inf (1-q)^(1-s) g(s;x) = (q^s;q)_inf * sum_n x^n / prod_{j<n} [s+j]_q,
# with the limit form at the pole s = 0: h(0;x) = (q;q)_inf (1-q) sum_{n>=1} x^n / prod_{1<=j<n} [j]_q.
# The factors (1-q)^(1-s) cancel between g(mu+a)g(mu+b) and g(mu)g(mu+a+b) because the shifts add up equally,
# and (q;q)_inf^2 > 0 is common, so sign(Delta) = sign(h(-2/5) h(0) - h(-7/10) h(3/10)).
# Every quantity is an interval [lo, hi] of rationals; products use all four corners and are rounded outward.
from fractions import Fraction as F
import sys
D = 1 << 160
def dn(v): return F(v.numerator * D // v.denominator, D)
def up(v): return F(-((-v.numerator * D) // v.denominator), D)
def imul(a, b):
    c = [a[0]*b[0], a[0]*b[1], a[1]*b[0], a[1]*b[1]]
    return (dn(min(c)), up(max(c)))
def iadd(a, b): return (a[0] + b[0], a[1] + b[1])
def pt(v): return (v, v)
r = F(93, 100); q = r**10
def qpow(t): return r**t                   # q^(t/10), exact
def br(t): return (1 - qpow(t)) / (1 - q)  # [t/10]_q, exact
def poch(t, N=120):
    # (q^(t/10); q)_inf for t > 0: product of the first N factors (each in (0,1)) times the tail factor
    # prod_{j>=N} (1 - a q^j) in [1 - a q^N/(1-q), 1]
    a = qpow(t); P = pt(F(1)); x = a
    for _ in range(N):
        P = imul(P, pt(1 - x)); x = x * q
    return imul(P, (1 - x / (1 - q), F(1)))
def pref(t):
    # (q^(t/10); q)_inf for t not in -10N; negative t: (1 - q^(t/10)) * (q^(t/10 + 1); q)_inf
    if t > 0: return poch(t)
    return imul(pt(1 - qpow(t)), poch(t + 10))
def series(t, X, K=160):
    # sum_{n>=0} c_n, c_0 = 1, c_{n+1} = c_n * X / [(t+10n)/10]_q.  Term intervals are propagated exactly
    # (outward rounding only).  From n = K on all terms share the sign of c_K and |c_{n+1}/c_n| <= rho
    # with rho = X / [(t+10K)/10]_q < 1 ([.]_q increasing), so the tail lies between 0 and c_K/(1-rho).
    S = pt(F(0)); c = pt(F(1))
    for n in range(K):
        S = iadd(S, c); c = imul(c, pt(X / br(t + 10 * n)))
    rho = X / br(t + 10 * K)
    assert 0 < rho < 1 and t + 10 * K > 0
    tail = imul(c, pt(1 / (1 - rho)))
    return (S[0] + min(F(0), tail[0]), S[1] + max(F(0), tail[1]))
X = F(int(sys.argv[1]) if len(sys.argv) > 1 else 93, 100)
assert 0 < X < 1 / (1 - q)                 # strictly inside the radius 1/(1-q)
hA = imul(pref(-4), series(-4, X))         # s = mu + alpha = -2/5
hB = imul(imul(poch(10), series(10, X)), pt(X * (1 - q)))   # s = mu + beta = 0 (limit form)
hC = imul(pref(-7), series(-7, X))         # s = mu = -7/10
hD = imul(pref(3), series(3, X))           # s = mu + alpha + beta = 3/10
A = imul(hA, hB); B = imul(hC, hD)
print('X', X, ' h(-2/5)h(0) in [%.9f, %.9f]   h(-7/10)h(3/10) in [%.9f, %.9f]' % (float(A[0]), float(A[1]), float(B[0]), float(B[1])))
print('certified Delta < 0:', A[1] < B[0])
