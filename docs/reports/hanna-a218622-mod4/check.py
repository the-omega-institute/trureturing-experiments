#!/usr/bin/env python3
"""Exact, reproducible checks for Hanna's A183161/A218622 conjecture.

OEIS data: Paul D. Hanna / OEIS Foundation, CC BY-SA 4.0; see README.md.
Python code: exact integers, Fraction and SymPy GF(2) polynomial arithmetic.
Finite checks are evidence, not a Lean proof of the conjecture.
"""
from collections import Counter
from fractions import Fraction
from math import comb
from pathlib import Path
import hashlib
import sympy as sp

x, y, g, k = sp.symbols('x y g k')


def poly(expression, modulus):
    terms = sp.Poly(sp.expand(expression), x, y).terms()
    return {ij: int(c) % modulus for ij, c in terms if int(c) % modulus}


def add(a, b, modulus):
    c = dict(a)
    for ij, value in b.items():
        c[ij] = (c.get(ij, 0) + value) % modulus
    return {ij: value for ij, value in c.items() if value}


def multiply(a, b, modulus):
    c = {}
    for (i, j), v in a.items():
        for (u, w), z in b.items():
            ij = i + u, j + w
            c[ij] = (c.get(ij, 0) + v * z) % modulus
    return {ij: value for ij, value in c.items() if value}


def section(a, r):
    return {(i // 2, j // 2): v for (i, j), v in a.items()
            if i % 2 == r and j % 2 == r}


def dense_transition(a, denominator, r, modulus):
    # Independent dense integer convolution: no sparse multiply or section.
    im = max((ij[0] for ij in a), default=0)
    jm = max((ij[1] for ij in a), default=0)
    iu = max(ij[0] for ij in denominator)
    ju = max(ij[1] for ij in denominator)
    left = [[a.get((i, j), 0) for j in range(jm + 1)] for i in range(im + 1)]
    right = [[denominator.get((i, j), 0) for j in range(ju + 1)]
             for i in range(iu + 1)]
    product = [[0] * (jm + ju + 1) for _ in range(im + iu + 1)]
    for i in range(im + 1):
        for j in range(jm + 1):
            for u in range(iu + 1):
                for v in range(ju + 1):
                    product[i + u][j + v] += left[i][j] * right[u][v]
    answer = {}
    for i in range(r, im + iu + 1, 2):
        for j in range(r, jm + ju + 1, 2):
            value = product[i][j] % modulus
            if value:
                answer[(i // 2, j // 2)] = value
    return answer


# The proposed certificate is an input to the checker, not its generated output.
ROWS = [
 (1+y+x*(2+y+y**3)+x**2*(1+2*y+2*y**2+3*y**3+y**4+3*y**7),
  y+x*y**2+x**2*(y**3+y**5), 1, 1, 2),
 (1+x*(2+y+y**3)+x**2*(1+y+2*y**3+y**4+y**5+2*y**6),
  y+x*(y+y**3)+x**2*(y**2+y**5)+x**3*(y**6+y**7), 1, 1, 3),
 (1+y+x*(3+2*y**2+3*y**4), y+x*(y+y**2+y**3)+x**2*y**4, 1, 4, 5),
 (3+3*y+x*(1+2*y**2+y**4), 1+y+x**2*y**3+x**3*y**7, 1, 6, 5),
 (1+x*(3+3*y+y**2+y**3), y+x*y**2+x**2*(y**3+y**5), 1, 7, 8),
 (x*(2+2*y**4), 1+x**2*y**5, 2, 9, 10),
 (3+x*(1+y+3*y**2+3*y**3),
  1+y+x*y**2+x**2*(y**2+y**3+y**4+y**5)+x**3*y**6, 1, 11, 12),
 (1+x*(3+y+y**2+3*y**3),
  y+x*(y+y**3)+x**2*(y**2+y**5)+x**3*(y**6+y**7), 1, 7, 12),
 (1+3*y+x*(1+2*y**2+y**4), y+x*(y+y**2+y**3)+x**2*y**4, 1, 4, 5),
 (0, 1+x*(y+y**2+y**3)+x**2*(y**3+y**4)+x**3*y**7, 2, 13, 14),
 (0, 1+x**2*y**5, 2, 9, 10),
 (3+x*(1+3*y+3*y**2+y**3), 1+y+x*(y+y**3)+x**2*(y**4+y**5)+x**3*y**7, 1, 11, 8),
 (3+y+x*(3+2*y**2+3*y**4), 1+y+x**2*y**3+x**3*y**7, 1, 6, 5),
 (0, 1+x**2*(y**2+y**4)+x**3*y**6, 2, 13, 9),
 (0, 0, 0, 14, 14),
]
Q = 1-x*(1+y)**3
D = 1+y+x*y+x**2*y+x**2*y**2+x**2*y**3+x**2*y**4+x**3*y**6+x**4*y**9
R0 = y+x*y**2+x**2*y**3+x**2*y**5
P0 = (Q+y-3*x*y**2*(1+y)**2)*Q
STATES = [(poly(a, 4), poly(b, 2)) for a, b, *_ in ROWS]
assert STATES[0] == (poly(P0, 4), poly(R0, 2))
Q2, D2 = poly(Q**2, 4), poly(D, 2)
assert poly(Q**4-Q.subs({x: x**2, y: y**2}, simultaneous=True)**2, 4) == {}
assert poly(D**2-D.subs({x: x**2, y: y**2}, simultaneous=True), 2) == {}
for index, (a, b, out, next0, next1) in enumerate(ROWS):
    ap, bp = STATES[index]
    assert (ap.get((0, 0), 0) + 2*bp.get((0, 0), 0)) % 4 == out
    assert out in (0, 1, 2)
    for r, nxt in enumerate((next0, next1)):
        sparse = section(multiply(ap, Q2, 4), r), section(multiply(bp, D2, 2), r)
        dense = dense_transition(ap, Q2, r, 4), dense_transition(bp, D2, r, 2)
        assert sparse == dense == STATES[nxt], (index, r)
    print(f'STATE {index:2}: output={out}; next0={next0}; next1={next1}; sparse=dense=OK')
print('CERTIFICATE: 15 input pairs, all 30 polynomial transition identities verified twice')

# Verify the elimination and rational-diagonal numerators, without numerical roots.
E = x**4*k**6+x**3*k**4+x**2*k**3+(x**2+1)*k**2+(x**2+x+1)*k+x**2
res = sp.resultant(x*g**3+g+1, g**4+g**2+g*k+k**2, g)
assert sp.Poly(res-E, x, k, modulus=2).is_zero
assert poly(sp.expand(E.subs({x: x*y, k: y}, simultaneous=True)/y)-D, 2) == {}
assert poly(y*sp.diff(E, k).subs({x: x*y, k: y}, simultaneous=True)-R0, 2) == {}
assert sp.Poly(sp.diff(E, k), x, k, modulus=2).eval({x: 0, k: 0}) == 1
print('ALGEBRA: exact GF(2) resultant, derivative, and denominator identities OK')


def automaton(n):
    state = 0
    while n:
        state = ROWS[state][3 + n % 2]
        n //= 2
    return ROWS[state][2]


def read_data(name):
    answer = {}
    for line in (Path(__file__).parent / name).read_text().splitlines():
        if line.strip() and not line.lstrip().startswith('#'):
            i, v = map(int, line.split())
            answer[i] = v
    return answer


published = read_data('b183161.txt')
residues = read_data('b218622.txt')
assert set(published) == set(range(301))
assert set(residues) == set(range(2049))
N = 300
convolution_root = [1]
source_sum = [1]
trinomial = [sum(comb(n, 2*j)*comb(2*j, j) for j in range(n//2+1)) for n in range(N+1)]
for n in range(1, N+1):
    s = sum(comb(n+j, n-j)*comb(2*n-j, j) for j in range(n+1))
    twice = s-sum(convolution_root[j]*convolution_root[n-j] for j in range(1, n))
    assert twice % 2 == 0
    convolution_root.append(twice//2)
    source_sum.append(sum((Fraction(trinomial[j]*comb(3*n-j, n-j)*2*j, 3*n-j)
                           for j in range(n+1)), Fraction(0)))
for n in range(N+1):
    assert source_sum[n].denominator == 1 if isinstance(source_sum[n], Fraction) else True
    assert published[n] == source_sum[n] == convolution_root[n]
    assert published[n] % 4 == residues[n] == automaton(n)
assert all(automaton(n) == v for n, v in residues.items())
print('FIDELITY: literal rational sum = normalized convolution = OEIS exact a(n), n=0..300')
print('FIDELITY: automaton = OEIS A218622 b-file residues, n=0..2048')
print('PREFIX:', ','.join(map(str, convolution_root[:12])))


def series_mul(a, b, limit, modulus=None):
    answer = [0]*(limit+1)
    for n in range(limit+1):
        answer[n] = sum(a[i]*b[n-i] for i in range(n+1))
        if modulus:
            answer[n] %= modulus
    return answer


# Integral source-series bridge and single binary correction, on a finite prefix.
M = 150
G = []
for n in range(M+1):
    value, remainder = divmod(comb(3*n, n), 2*n+1)
    assert remainder == 0
    G.append(value)
t = G.copy()
t[0] = 0
t2 = series_mul(t, t, M)
power = [1]+[0]*M
central_sum = [0]*(M+1)
for m in range(M//2+1):
    central = comb(2*m, m)
    central_sum = [a+central*b for a, b in zip(central_sum, power)]
    power = series_mul(power, t2, M)
F = series_mul(G, central_sum, M)
assert F == convolution_root[:M+1]
g2 = series_mul(G, G, M)
g4 = series_mul(g2, g2, M)
radicand = [int(n == 0) - (2*g2[n-1] if n else 0) - (3*g4[n-2] if n >= 2 else 0)
            for n in range(M+1)]
assert series_mul(series_mul(F, F, M), radicand, M) == [1]+[0]*M
h = [0]*(M+1)
p = 2
while p <= M:
    for n in range(p, M+1, p):
        h[n] ^= G[n//p] % 2
    p *= 2
correction = series_mul([v % 2 for v in G], h, M, 2)
assert all((G[n]+2*correction[n]) % 4 == F[n] % 4 for n in range(M+1))
for m in range(1, 1001):
    assert comb(2*m, m) % 4 == (2 if m & (m-1) == 0 else 0)
print('BRIDGE: F=G*Sum binom(2m,m)*(G-1)^(2m); source square-root identity and G+2k, n=0..150')
print('CENTRAL BINOMIAL: predicted mod-4 power-of-two support, m=1..1000')


def rational_diagonal(numerator, denominator, modulus, limit):
    assert denominator.get((0, 0)) == 1
    terms = [(i, j, value) for (i, j), value in denominator.items() if (i, j) != (0, 0)]
    coefficients = [[0]*(limit+1) for _ in range(limit+1)]
    for i in range(limit+1):
        for j in range(limit+1):
            coefficients[i][j] = (numerator.get((i, j), 0)-sum(
                value*coefficients[i-u][j-v] for u, v, value in terms if u <= i and v <= j)) % modulus
    return [coefficients[i][i] for i in range(limit+1)]


diagG = rational_diagonal(STATES[0][0], Q2, 4, 64)
diagK = rational_diagonal(STATES[0][1], D2, 2, 64)
assert diagG == [v % 4 for v in G[:65]]
assert diagK == correction[:65]
assert [(a+2*b) % 4 for a, b in zip(diagG, diagK)] == [automaton(n) for n in range(65)]
print('DIAGONALS: independent bivariate rational-series division agrees, n=0..64')
print('AUTOMATON COUNTS n=0..65535:', dict(sorted(Counter(automaton(n) for n in range(65536)).items())))
print('SCOPE: finite source checks and finite polynomial certificates; no Lean proof executed')
print('check.py sha256:', hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
print('ALL_CHECKS_OK')
