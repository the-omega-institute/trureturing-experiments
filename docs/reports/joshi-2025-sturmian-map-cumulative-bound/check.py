"""Check the "stronger cumulative bound on A(d)" of G. Joshi's PhD thesis (Open University 2025,
DOI 10.21954/ou.ro.00106631), Section 6.1, p. 118:

    for every n in N,  d <= q_n  ==>  A(d) <= A(q_n)   "could be true for all g(d)"

for the Sturmian sequence s of slope 0 < alpha < 1/2 (Definition 2.2.3: s_n = b iff
<(n+1)alpha> lies in [1 - alpha, 1)), where A(d) is the longest monochromatic arithmetic
progression of difference d over both letters and every start (Definition 2.6.2) and q_n is the
n-th convergent denominator of alpha.  Proposition 4.4.4 of the thesis proves it when
g(d) = min(<d alpha>, 1 - <d alpha>) <= alpha.

Prediction registered before the run:
  (1) alpha = 2 - sqrt(3) = [0; 3, 1, 2, 1, 2, ...] violates it at n = 1: q_1 = 3, A(2) = 8 > A(3) = 4
      (exact integer certificates for the witness MAP and for the bound A(3) <= 4);
  (2) the thesis formulae of Corollary 4.3.13 reproduce the brute-force values (and the thesis
      examples A(5) = 5, A(3) = 6 for <sqrt 5>, A(3) = 5 for the Fibonacci word);
  (3) no violation occurs among differences with g(d) <= alpha (Proposition 4.4.4) in any scan;
  (4) every violation found has g(d) > alpha (the Case-2 cluster regime of Theorem 4.3.9).
Negative controls: (2) and (3), and the thesis's own running slopes <sqrt 5> and the Fibonacci
slope, for which no violation is expected in the scanned range.
Usage: python3 check.py   (standard library only; about a minute)
"""
import sys
from math import isqrt, ceil, floor
from decimal import Decimal, getcontext

getcontext().prec = 60
OK = True


def fail(msg):
    global OK
    OK = False
    print("FAIL", msg)


class Quad:
    """alpha = (P + Q*sqrt(D)) / R with integers, D not a square, R > 0."""

    def __init__(self, P, Q, D, R):
        self.P, self.Q, self.D, self.R = P, Q, D, R

    def floor_mul(self, n):
        t = n * self.Q
        s = isqrt(t * t * self.D) if t >= 0 else -(isqrt(t * t * self.D) + 1)
        return (n * self.P + s) // self.R

    def dec(self):
        return (Decimal(self.P) + Decimal(self.Q) * Decimal(self.D).sqrt()) / Decimal(self.R)


def word_quad(a, N):
    return [a.floor_mul(n + 2) - a.floor_mul(n + 1) for n in range(N)]


def word_dec(x, N):
    return [int(floor((n + 2) * x)) - int(floor((n + 1) * x)) for n in range(N)]


def A_brute(w, d):
    best = 0
    for r in range(d):
        run, prev = 0, None
        for i in range(r, len(w), d):
            if w[i] == prev:
                run += 1
            else:
                run, prev = 1, w[i]
            if run > best:
                best = run
    return best


def cf_dens(x, k):
    a, y = [], x
    for _ in range(k):
        y = 1 / y
        ai = int(floor(y))
        a.append(ai)
        y = y - ai
    q, qm = [1], 0
    for ai in a:
        q.append(ai * q[-1] + qm)
        qm = q[-2]
    return a, q


def frac(x):
    return x - int(floor(x))


def g_of(x, d):
    f = frac(d * x)
    return min(f, 1 - f)


def A_formula(x, d):
    """Corollary 4.3.13 of the thesis (0 < alpha < 1/2)."""
    g = g_of(x, d)
    if g <= x:
        return int(ceil((1 - x) / g))
    k = int(ceil((1 - x) / g))
    if g < Decimal(1) / k:
        return k * int(ceil((1 - x - (k - 1) * g) / (1 - k * g)))
    return k * int(ceil((g - x) / (k * g - 1))) + k - 2


# ---------------------------------------------------------------- (1) the exact witness
print("== (1) alpha = 2 - sqrt(3)")
a0 = Quad(2, -1, 3, 1)
# 0 < alpha < 1/2: 1 < sqrt3 < 2 and sqrt3 > 3/2 (9 > 4*... ) checked by integer squares
if not (3 > 1 and 3 < 4 and 4 * 3 > 9):
    fail("range of alpha")
# q_1 = floor(1/alpha) = floor(2 + sqrt 3) = 3 since 1 < sqrt 3 < 2
q1 = 2 + isqrt(3)
print("q_1 =", q1)
if q1 != 3:
    fail("q_1")
w0 = word_quad(a0, 4000)
# witness: positions 1,3,5,7,9 (difference 2) all letter a (= 0); <m alpha> < 1 - alpha for m = 2..10 even
wit = [w0[1 + 2 * k] for k in range(5)]
print("letters at 1,3,5,7,9:", wit)
if wit != [0] * 5:
    fail("witness letters")
# the same five facts as exact integer inequalities: <m alpha> = ceil(m sqrt3) - m sqrt3 < sqrt3 - 1
#   <=>  (m+1) sqrt3 > ceil(m sqrt3) + 1  <=>  3 (m+1)^2 > (isqrt(3 m^2) + 2)^2
for m in (2, 4, 6, 8, 10):
    if not (3 * (m + 1) ** 2 > (isqrt(3 * m * m) + 2) ** 2):
        fail(f"exact inequality m={m}")
# A(3) <= 4: g = 1 - 3 alpha = 3 sqrt3 - 5; need 0 < g < alpha, 4 g > 1 - alpha, 2 g > alpha
checks = {
    "g > 0      (27 > 25)": 27 > 25,
    "g < alpha  (48 < 49)": 48 < 49,
    "4g > 1-alpha (363 > 361)": 363 > 361,
    "2g > alpha (147 > 144)": 147 > 144,
}
for k_, v in checks.items():
    print(" ", k_, v)
    if not v:
        fail(k_)
x0 = a0.dec()
for k_, v in {"g>0": 3 * Decimal(3).sqrt() - 5 > 0, "g<alpha": 3 * Decimal(3).sqrt() - 5 < x0,
              "4g>1-alpha": 4 * (3 * Decimal(3).sqrt() - 5) > 1 - x0,
              "2g>alpha": 2 * (3 * Decimal(3).sqrt() - 5) > x0}.items():
    if not v:
        fail("decimal cross-check " + k_)
A2, A3 = A_brute(w0, 2), A_brute(w0, 3)
print(f"brute force on 4000 letters: A(1)={A_brute(w0,1)} A(2)={A2} A(3)={A3}")
if (A2, A3) != (8, 4):
    fail("brute A(2), A(3)")
first8 = next(i for i in range(len(w0) - 16) if all(w0[i + 2 * k] == 0 for k in range(8)))
print("first 8-term a-MAP of difference 2 starts at", first8)

# ---------------------------------------------------------------- (2) formula vs brute force
print("== (2) Corollary 4.3.13 vs brute force")
D5 = Decimal(5).sqrt()
samples = {
    "<sqrt5>": Quad(-2, 1, 5, 1),
    "fibonacci (3-sqrt5)/2": Quad(3, -1, 5, 2),
    "2-sqrt3": a0,
    "(14-sqrt26)/34": Quad(14, -1, 26, 34),
}
words = {nm: word_quad(a, 120000) for nm, a in samples.items()}
for nm, a in samples.items():
    x = a.dec()
    mism = [(d, A_formula(x, d), A_brute(words[nm], d)) for d in range(1, 41) if A_formula(x, d) != A_brute(words[nm], d)]
    print(f"  {nm}: d=1..40 mismatches {mism}")
    if mism:
        fail("formula mismatch " + nm)
xs = samples["<sqrt5>"].dec()
if (A_formula(xs, 5), A_formula(xs, 3)) != (5, 6):
    fail("thesis Examples 4.3.6/4.3.10")
if A_brute(words["fibonacci (3-sqrt5)/2"], 3) != 5:
    fail("thesis Example 2.6.5")
print("  thesis examples: <sqrt5> A(5)=5, A(3)=6; Fibonacci A(3)=5: reproduced")

# ---------------------------------------------------------------- (3),(4) violation scans
print("== (3)(4) violations of d <= q_n => A(d) <= A(q_n), n >= 1, q_n <= 250")


def scan(x, qmax=250):
    a, q = cf_dens(x, 30)
    out = []
    for n in range(1, len(q)):
        if q[n] > qmax:
            break
        Aq = A_formula(x, q[n])
        for d in range(1, q[n] + 1):
            Ad = A_formula(x, d)
            if Ad > Aq:
                out.append((n, q[n], d, Ad, Aq, g_of(x, d) > x))
    return a, q, out


for nm, a in samples.items():
    x = a.dec()
    cf, q, v = scan(x)
    lv = sorted({t[0] for t in v})
    print(f"  {nm}: cf {cf[:8]} q {q[:7]} violations {len(v)} at levels n={lv}")
    for t in v[:4]:
        print("     n=%d q_n=%d d=%d A(d)=%d A(q_n)=%d g(d)>alpha=%s" % t)
    if any(not t[5] for t in v):
        fail("violation with g(d) <= alpha (contradicts Proposition 4.4.4) " + nm)
    if nm in ("<sqrt5>", "fibonacci (3-sqrt5)/2") and v:
        fail("unexpected violation for a thesis running slope " + nm)
    if nm == "2-sqrt3" and (1, 3, 2, 8, 4, True) not in v:
        fail("expected witness missing from scan")
# brute-force confirmation of the formula-detected violations for (14 - sqrt 26)/34
x = samples["(14-sqrt26)/34"].dec()
_, _, v = scan(x)
wlong = word_quad(samples["(14-sqrt26)/34"], 1200000)
conf = 0
for t in v:
    n, qn, d, Ad, Aq, _ = t
    if 3 * max(Ad * d, Aq * qn) > len(wlong):
        continue  # the prefix is too short to exhibit the maximum
    conf += 1
    if (A_brute(wlong, d), A_brute(wlong, qn)) != (Ad, Aq):
        fail(f"brute confirmation {t}")
print(f"  brute force on 1.2e6 letters confirms {conf} of the {len(v)} violations of (14-sqrt26)/34 (the rest need a longer prefix)")
if conf < 4:
    fail("too few brute confirmations")

# family alpha_B = [0; 3, (1, B, 1)] : ratio max_{d<=q_n} A(d) / A(q_n) along the levels
print("== family [0;3,(1,B,1)]: worst ratio max_{d<=q_n}A(d)/A(q_n) per level (q_n <= 2000)")


def cf_value(prefix, period, reps=40):
    seq = prefix + period * reps
    y = Decimal(0)
    for t in reversed(seq):
        y = 1 / (t + y)
    return y


for B in (4, 6, 8, 12):
    x = cf_value([3], [1, B, 1])
    a, q = cf_dens(x, 12)
    rows = []
    for n in range(1, len(q)):
        if q[n] > 2000:
            break
        Aq = A_formula(x, q[n])
        best = max(range(1, q[n] + 1), key=lambda d: A_formula(x, d))
        rows.append((n, q[n], best, A_formula(x, best), Aq, round(A_formula(x, best) / Aq, 3)))
    print(f"  B={B}:", [(r[0], r[1], r[2], r[3], r[4], r[5]) for r in rows if r[5] > 1])

# random slopes: how often does level n = 1..3 fail?
print("== random slopes alpha uniform in (0, 1/2), levels with q_n <= 60")
import random
random.seed(20261016)
cnt = viol = 0
for _ in range(400):
    x = Decimal(random.random()) / 2
    if x <= Decimal("0.001"):
        continue
    _, q, v = scan(x, qmax=60)
    cnt += 1
    viol += bool(v)
print(f"  {viol}/{cnt} random slopes violate at some level with q_n <= 60")
# level n = 1 on (1/4, 1/3): A(2) = 2 ceil(alpha/(4 alpha - 1)), A(3) = ceil((1-alpha)/(1-3 alpha))
grid = [Decimal(1) / 4 + Decimal(i) / 100000 * (Decimal(1) / 3 - Decimal(1) / 4) for i in range(1, 100000)]
bad = [x for x in grid if A_formula(x, 2) > A_formula(x, 3)]
print(f"  level 1 on (1/4,1/3): violating fraction of a 1e5 grid = {len(bad)/len(grid):.4f}; sup of violating grid points = {float(max(bad)):.6f}")
print("ALL_OK" if OK else "SOME_FAILED")
sys.exit(0 if OK else 1)
