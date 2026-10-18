"""OEIS A163921 (P. D. Hanna) conjecture and OEIS A157675.

Claim: if F(x) = exp(sum_{n>=1} L(n) x^n/n) has integer coefficients, then the solution G (G(0) = 1) of
    G(x) = exp( sum_{n>=1} G(L(n) x)^n x^n / n )
has integer coefficients. A163921 is L = sigma, A157675 is L(n) = 2^n.
"""
import sys, random
from fractions import Fraction as Fr

N = int(sys.argv[1]) if len(sys.argv) > 1 else 16
ok = True
def check(name, cond):
    global ok
    print(("PASS " if cond else "FAIL ") + name); ok = ok and bool(cond)

def mul(a, b):
    c = [Fr(0)] * N
    for i, x in enumerate(a):
        if x:
            for j in range(N - i):
                c[i + j] += x * b[j]
    return c

def exp_series(h):            # h[0] == 0; E' = h' E
    e = [Fr(0)] * N; e[0] = Fr(1)
    for n in range(1, N):
        e[n] = sum(k * h[k] * e[n - k] for k in range(1, n + 1)) / n
    return e

def powr(a, m):
    r = [Fr(0)] * N; r[0] = Fr(1)
    for _ in range(m):
        r = mul(r, a)
    return r

def solve(L):
    G = [Fr(0)] * N; G[0] = Fr(1)
    for _ in range(N):          # coefficient n of the output depends only on coefficients < n of the input
        H = [Fr(0)] * N
        for n in range(1, N):
            Gs = [G[k] * Fr(L(n)) ** k for k in range(N)]   # G(L(n) x)
            P = powr(Gs, n)
            for k in range(N - n):
                H[n + k] += P[k] / n
        G = exp_series(H)
    return G

def F_of(L):
    return exp_series([Fr(0)] + [Fr(L(n), n) for n in range(1, N)])

integral = lambda s: all(x.denominator == 1 for x in s)
sigma = lambda n: sum(d for d in range(1, n + 1) if n % d == 0)

G = solve(sigma)
A163921 = [1, 1, 2, 7, 41, 385, 5769]
check(f"L = sigma: F integral {integral(F_of(sigma))}; G integral to x^{N-1}; G starts {[int(x) for x in G[:7]]}",
      integral(F_of(sigma)) and integral(G) and [int(x) for x in G[:7]] == A163921)
G2 = solve(lambda n: 2 ** n)
A157675 = [1, 1, 3, 19, 237, 5741, 270857, 25099497]
check(f"L = 2^n: F = 1/(1-2x) integral; G integral to x^{N-1}; G starts {[int(x) for x in G2[:8]]}",
      integral(F_of(lambda n: 2 ** n)) and integral(G2) and [int(x) for x in G2[:8]] == A157675)

rng = random.Random(20261018)
for t in range(4):
    c = {d: rng.randint(-3, 3) for d in range(1, N)}         # F = prod (1 - x^d)^(-c_d)
    L = lambda n, c=c: sum(d * c[d] for d in range(1, n + 1) if n % d == 0)
    Fi = integral(F_of(L)); Gi = integral(solve(L))
    check(f"random integral F #{t} (c_1..c_5 = {[c[d] for d in range(1, 6)]}): F integral {Fi}, G integral {Gi}", Fi and Gi)

Lbad = lambda n: 1 if n == 2 else 0                            # F = exp(x^2/2), not integral
Gb = solve(Lbad)
check(f"control L = [n = 2]: F not integral, G not integral (first non-integer at x^{next(i for i, x in enumerate(Gb) if x.denominator != 1)})",
      not integral(F_of(Lbad)) and not integral(Gb))
print("RESULT", "PASS" if ok else "FAIL"); sys.exit(0 if ok else 1)
