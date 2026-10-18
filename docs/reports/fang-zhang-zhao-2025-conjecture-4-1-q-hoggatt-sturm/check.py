"""Fang-Zhang-Zhao, JMAA 561 (2026) 130615 (arXiv:2505.05873v1), Conjecture 4.1: for every integer q > 0 and
m = 2, 3 the q,t-Hoggatt polynomials H_n^[m](q, t) = sum_{k+l=n-1} M(k,l,m;q) t^k form a Sturm sequence
(consecutive members strictly interlace, all roots negative), where
M(k,l,m;q) = q^{m k(k+1)/2} * sum over plane partitions in a k x l x m box of q^{|pi|}."""
import sys, itertools
from fractions import Fraction as Fr
import sympy as sp

ok = True
def check(name, cond):
    global ok
    print(("PASS " if cond else "FAIL ") + name); ok = ok and bool(cond)

def pp_gf(k, l, m, q):            # brute force: k x l arrays, entries 0..m, weakly decreasing along rows and columns
    if k == 0 or l == 0: return Fr(1)
    total = Fr(0)
    rows = [r for r in itertools.product(range(m, -1, -1), repeat=l) if all(r[i] >= r[i+1] for i in range(l-1))]
    def rec(i, prev, s):
        nonlocal total
        if i == k: total += Fr(q) ** s; return
        for r in rows:
            if all(r[j] <= prev[j] for j in range(l)): rec(i + 1, r, s + sum(r))
    rec(0, (m,) * l, 0); return total

def qint(r, q): return sum(Fr(q) ** j for j in range(r))
def qfact(r, q):
    p = Fr(1)
    for j in range(1, r + 1): p *= qint(j, q)
    return p
def qbin(N, k, q): return qfact(N, q) / (qfact(k, q) * qfact(N - k, q)) if 0 <= k <= N else Fr(0)
def macmahon(k, l, m, q):
    p = Fr(1)
    for r in range(m):
        p *= qfact(k + l + r, q) * qfact(r, q) / (qfact(k + r, q) * qfact(l + r, q))
    return p

good = all(pp_gf(k, l, m, q) == macmahon(k, l, m, q) for q in (2, 3) for m in (2, 3) for k in range(0, 4) for l in range(0, 4))
check("MacMahon box product equals the brute-force plane-partition sum (q = 2, 3; m = 2, 3; k, l <= 3)", good)

def H(n, m, q):                   # coefficients of H_n^[m] in t, via MacMahon
    return [Fr(q) ** (m * k * (k + 1) // 2) * macmahon(k, n - 1 - k, m, q) for k in range(n)]
def F(N, m, q): return [Fr(q) ** (m * k * (k + 1) // 2) * prod_qbin(N, k, m, q) for k in range(N - m + 2)]
def prod_qbin(N, k, m, q):
    p = Fr(1)
    for i in range(m): p *= qbin(N, k + i, q)
    return p
def prod_j(N, m, q):
    p = Fr(1)
    for j in range(1, m): p *= qbin(N, j, q)
    return p
lem1 = all(H(n, m, q) == [c / prod_j(n + m - 2, m, q) for c in F(n + m - 2, m, q)] for q in (2, 3) for m in (2, 3) for n in range(1, 8))
check("Lemma 1: H_n^[m] = F_{n+m-2,m} / prod_{j<m} C_q(n+m-2, j) (q = 2, 3; m = 2, 3; n <= 7)", lem1)

def had(a, b): return [x * y for x, y in zip(a, b)]
def Dq(a, q): return [qint(k, q) * a[k] for k in range(1, len(a))]
def lem5(N, m, q):
    G = had(F(N, 1, q), F(N - 1, m - 1, q)); lhs = Dq(G, q)
    rhs = [qint(N, q) * Fr(q) ** m * c * Fr(q) ** (m * k) for k, c in enumerate(F(N - 1, m, q))]
    return lhs == rhs[:len(lhs)] and len(lhs) == len(rhs)
check("Lemma 5: D_q(F_{N,1} * F_{N-1,m-1})(t) = [N]_q q^m F_{N-1,m}(q^m t) (q = 2, 3; m = 2, 3, 4; N <= 9)",
      all(lem5(N, m, q) for q in (2, 3) for m in (2, 3, 4) for N in range(m, 10)))

t = sp.symbols('t')
def roots(coeffs):
    P = sp.Poly(list(reversed([sp.Rational(c.numerator, c.denominator) for c in coeffs])), t)
    rr = sp.real_roots(P)
    return rr if len(rr) == P.degree() else None
def strict_interlace(f, g):       # deg g = deg f + 1, all roots negative, y1 < x1 < y2 < ... < x_d < y_{d+1}
    x, y = roots(f), roots(g)
    if x is None or y is None: return False
    if any(sp.N(v) >= 0 for v in list(x) + list(y)): return False
    seq = []
    for i in range(len(x)): seq += [y[i], x[i]]
    seq.append(y[-1])
    return all(sp.N(seq[i] - seq[i + 1], 50) < 0 for i in range(len(seq) - 1))
res = {}
for q in (Fr(3, 2), Fr(2), Fr(3)):
    for m in (2, 3, 4):
        res[(q, m)] = all(strict_interlace(H(n, m, q), H(n + 1, m, q)) for n in range(1, 11))
check("strict consecutive interlacing H_n -> H_{n+1}, n <= 10, for q in {3/2, 2, 3}, m in {2, 3, 4}: " + str({f'{k[0]},{k[1]}': v for k, v in res.items()}), all(res.values()))
print("RESULT", "PASS" if ok else "FAIL"); sys.exit(0 if ok else 1)
