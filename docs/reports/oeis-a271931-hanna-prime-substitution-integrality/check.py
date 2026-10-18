"""OEIS A271931/A271932/A271933 (Hanna 2016): for prime p there is an integer series G with
G(x)^p = G(x^p + p x G(x)^p), G(0) = 0, G'(0) = 1.

Route: the reversion R = G^(-1) satisfies R(y^p) = R(y)^p + p y^p R(y). Comparing coefficients of y^(n+p-1),
p r_n = [y^(n+p-1)] (R(y^p) - R_<n(y)^p - p y^p R_<n(y)), and the bracket is divisible by p since
R_<n(y)^p = R_<n(y^p) (mod p) and (n+p-1)/p < n. So R and G are integral.
Checks: exact divisions for primes p <= 23 (orders below); G matches the OEIS data for p = 5, 7, 11;
G^p = G(x^p + p x G^p) holds to the computed order; for composite p = 4, 6, 8, 9 the solution is not integral.
"""
from fractions import Fraction

def mul(a, b, n):
    c = [0] * (n + 1)
    for i, x in enumerate(a[:n + 1]):
        if x:
            for j, y in enumerate(b[:n + 1 - i]):
                if y: c[i + j] += x * y
    return c

def power(a, k, n):
    r = [1] + [0] * n
    for _ in range(k): r = mul(r, a, n)
    return r

def compose(f, g, n):
    res = [0] * (n + 1); P = [1] + [0] * n
    for k in range(1, n + 1):
        P = mul(P, g, n)
        if f[k]:
            for i in range(n + 1): res[i] += f[k] * P[i]
    return res

def reversion_coeffs(p, N, exact_div=True):
    zero = 0 if exact_div else Fraction(0)
    r = [zero] * (N + p + 1); r[1] = 1 if exact_div else Fraction(1); exact = True
    for n in range(2, N + 1):
        d = n + p - 1
        Rt = r[:n] + [zero] * (N + p + 1 - n)
        Rp = power(Rt, p, d)
        lhs = r[d // p] if d % p == 0 else zero
        pyR = sum(p * Rt[i] for i in range(1, n) if i + p == d)
        rest = lhs - Rp[d] - pyR
        if exact_div:
            exact &= rest % p == 0; r[n] = rest // p
        else:
            r[n] = Fraction(rest) / p
    return exact, r

def inverse(r, N):
    G = [0] * (N + 1); G[1] = 1
    for m in range(2, N + 1):
        G[m] = -compose(r, G, m)[m]
    return G

def oeis_terms(path, k):
    vals = []
    for line in open(path):
        if line[:2] in ('%S', '%T', '%U'):
            vals += [int(t) for t in line.split(None, 2)[2].strip().strip(',').split(',') if t]
    return vals[:k]

ok = True
orders = {2: 30, 3: 30, 5: 30, 7: 25, 11: 20, 13: 18, 17: 14, 19: 13, 23: 12}
for p, N in orders.items():
    exact, r = reversion_coeffs(p, N)
    G = inverse(r, N)
    Gp = power(G, p, N + p)
    inner = [0] * (N + p + 1); inner[p] = 1
    for i, v in enumerate(mul([0] + [p] + [0] * (N + p - 1), power(G, p, N + p), N + p)): inner[i] += v
    rhs = compose(G + [0] * p, inner, N + p - 1)
    eq = Gp[:N + 1] == rhs[:N + 1]
    print(f'p = {p:2d}: divisions exact {exact}, equation holds to order {N}: {eq}, G = {G[1:7]}')
    ok &= exact and eq
    if p in (5, 7, 11):
        k = min(N, 12); ref = oeis_terms(f'oeis-a2719{ {5: "31", 7: "32", 11: "33"}[p] }.txt', k)
        print(f'        matches OEIS A2719{ {5: "31", 7: "32", 11: "33"}[p] } (first {k} terms): {G[1:k + 1] == ref}')
        ok &= G[1:k + 1] == ref
for p in (4, 6, 8, 9):
    _, r = reversion_coeffs(p, 12, exact_div=False)
    bad = [n for n in range(1, 13) if r[n].denominator != 1]
    print(f'composite p = {p}: first non-integer reversion coefficient r_{bad[0]} = {r[bad[0]]}'); ok &= bool(bad)
print('RESULT', 'PASS' if ok else 'FAIL')
