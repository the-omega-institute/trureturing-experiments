#!/usr/bin/env python3
"""Exact checks for Hanna's OEIS A396804 modulo-five conjecture.

Python 3 standard library only. No floats, no Lean proof. The independent
integer Bell recurrence computes 201 terms; Fraction/Horner fixed-point
iteration cross-checks 17. The quotient-ring calculation checks a candidate
certificate, whose universal Hurwitz interpretation remains a proof obligation.
"""
from fractions import Fraction
from hashlib import sha256
from math import comb, factorial, prod
from urllib.request import Request, urlopen

N = 201
SOURCE_FIRST = [1, 2, 27, 820, 41005, 2933046, 277407235, 33008920904,
                4775699432745, 819537205428010, 163708079071865911,
                37510927241086166220, 9743252939007037322509,
                2841190409374953654166526, 922635217730641257310962075,
                331365807244074316193652673936,
                130847111194991759214337098786385]
BFILE_URL = 'https://oeis.org/A396804/b396804.txt'
BFILE_SHA256 = '32a0947573161ef08b5eb7b11c581d6c676bbcba1001a39a5634c599e55ac986'
TABLE = [0, 1, 0, 1, 0, 0, 1, 0, 4, 0, 0, 1, 0, 4, 1, 0, 1, 0, 0, 2]


def integer_fixed_point(limit):
    """Triangular Bell recursion, not iteration of the proposed certificate.

    g[k] is A composed k times. At stage n, a[n]=n*exp(g[4])[n-1]
    uses only earlier data. B_g(n,k)=sum_j C(n-1,j-1)g[j]B_g(n-j,k-1).
    """
    g = [[0] * (limit + 1) for _ in range(5)]
    bells = [[[0] * (limit + 1) for _ in range(limit + 1)] for _ in range(3)]
    for b in bells:
        b[0][0] = 1
    exp4 = [1] + [0] * limit
    for n in range(1, limit + 1):
        g[1][n] = n * exp4[n - 1]
        row = [comb(n - 1, j - 1) for j in range(1, n + 1)]
        for k in range(1, 4):
            b, inner = bells[k - 1], g[k]
            for r in range(1, n + 1):
                b[n][r] = sum(row[j - 1] * inner[j] * b[n - j][r - 1]
                              for j in range(1, n - r + 2))
            g[k + 1][n] = sum(g[1][r] * b[n][r] for r in range(1, n + 1))
        exp4[n] = sum(row[j - 1] * g[4][j] * exp4[n - j]
                      for j in range(1, n + 1))
    return g, exp4


def rational_fixed_point(limit):
    """Ordinary rational FPS via convolution, Horner composition and exp."""
    def mul(f, g):
        return [sum((f[i] * g[n - i] for i in range(n + 1)), Fraction(0))
                for n in range(limit + 1)]

    def compose(f, g):
        v = [Fraction(0)] * (limit + 1)
        for c in reversed(f):
            v = mul(v, g)
            v[0] += c
        return v

    f = [Fraction(0)] * (limit + 1)
    f[1] = 1
    for _ in range(limit):
        g = f
        for _ in range(3):
            g = compose(f, g)
        e = [Fraction(1)] + [Fraction(0)] * limit
        for n in range(1, limit + 1):
            e[n] = sum((j * g[j] * e[n - j] for j in range(1, n + 1)),
                       Fraction(0)) / n
        f = [Fraction(0)] + e[:-1]
    return [f[n] * factorial(n) for n in range(limit + 1)]


# R = F_5[X,E]/(X^5,E^5-1), represented in its 25-monomial basis.
def clean(p):
    return {k: v % 5 for k, v in p.items() if v % 5}


def add(*ps):
    out = {}
    for p in ps:
        for k, v in p.items():
            out[k] = out.get(k, 0) + v
    return clean(out)


def scale(c, p):
    return clean({k: c * v for k, v in p.items()})


def mul(p, q):
    out = {}
    for (r, j), a in p.items():
        for (s, k), b in q.items():
            if r + s < 5:
                key = (r + s, (j + k) % 5)
                out[key] = out.get(key, 0) + a * b
    return clean(out)


ONE, X, E = {(0, 0): 1}, {(1, 0): 1}, {(0, 1): 1}


def power(p, n):
    out = ONE
    for _ in range(n):
        out = mul(out, p)
    return out


def exp_r(g):
    assert all(r >= 1 for r, _ in g)
    assert sum(v for (r, _), v in g.items() if r == 1) % 5 == 1
    h = add(g, scale(-1, X))
    return mul(E, add(*(scale(c, power(h, r))
                        for r, c in enumerate([1, 1, 3, 1, 4]))))


def compose_r(p, g):
    eg = exp_r(g)
    return add(*(scale(c, mul(power(g, r), power(eg, j)))
                 for (r, j), c in p.items()))


def evaluate(p, n):
    # [gamma_n] X^r E^j = (n)_r * j^(n-r); 0^0=1.
    return sum(c * prod(range(n - r + 1, n + 1)) * pow(j, n - r, 5)
               for (r, j), c in p.items() if r <= n) % 5


def exp_hurwitz(g):
    out = [1] + [0] * (len(g) - 1)
    for n in range(1, len(g)):
        out[n] = sum(comb(n - 1, j - 1) * g[j] * out[n - j]
                     for j in range(1, n + 1)) % 5
    return out


def residue_formula(n):
    assert n >= 5
    fall = lambda r: prod(range(n - r + 1, n + 1))
    return (n - fall(2) + fall(3) * (-1 + pow(2, n - 3, 5) - pow(4, n - 3, 5))
            + fall(4) * (1 + pow(2, n - 4, 5) + 3 * pow(4, n - 4, 5))) % 5


def main():
    g, exp4 = integer_fixed_point(N)
    a = g[1]
    assert a[1:18] == SOURCE_FIRST
    assert rational_fixed_point(17) == a[:18]
    assert g[4][1:8] == [1, 8, 180, 7240, 435400, 35615124, 3738724696]
    print('Fidelity: integer Bell recursion and Fraction/Horner agree through n=17')
    print('Source first eight a(n):', a[1:9])
    print('Source fourth iterate through n=7:', g[4][1:8])
    data = urlopen(Request(BFILE_URL, headers={'User-Agent': 'A396804 exact verifier'}),
                   timeout=60).read()
    assert sha256(data).hexdigest() == BFILE_SHA256
    published = [tuple(map(int, s.split())) for s in data.decode('ascii').splitlines()
                 if s.strip() and not s.lstrip().startswith('#')]
    assert published == list(enumerate(a[1:], start=1))
    print('Fidelity: all 201 integer terms equal the pinned OEIS b-file')
    print('b-file sha256:', sha256(data).hexdigest())
    w6 = Fraction(g[4][6], factorial(5))
    assert w6 == Fraction(2967927, 10) and w6.denominator % 5 == 0
    print('Sixth logarithmic weight of A/x:', w6, '(not 5-integral)')

    f = {(1, 1): 1, (2, 0): 1, (2, 1): 4,
         (3, 0): 1, (3, 1): 4, (3, 2): 1, (3, 4): 4,
         (4, 1): 1, (4, 2): 1, (4, 4): 3}
    expected = [([1, 4, 0, 0, 0], [1, 4, 1, 0, 4]),
                ([0, 1, 4, 0, 0], [0, 2, 0, 2, 1]),
                ([0, 0, 1, 3, 1], [0, 1, 2, 2, 0]),
                ([1, 1, 0, 1, 2], [2, 2, 3, 0, 3])]
    assert power(X, 5) == {} and power(E, 5) == ONE
    h = X
    for k, (v, w) in enumerate(expected, 1):
        h = compose_r(f, h)
        assert [h.get((2, j), 0) for j in range(5)] == v
        assert [h.get((3, j), 0) for j in range(5)] == w
        vals = [evaluate(h, n) for n in range(N + 1)]
        assert vals == [c % 5 for c in g[k]]
        assert [evaluate(exp_r(h), n) for n in range(N + 1)] == exp_hurwitz(vals)
        print('Quotient recurrence k=%d: V=%s; W=%s' % (k, v, w))
    rhs = mul(X, exp_r(h))
    assert rhs == f
    print('Exact quotient identity: f = X Exp_R(f composed four times); all 25 coefficients agree')
    assert [evaluate(f, n) for n in range(N + 1)] == [c % 5 for c in a]
    assert [residue_formula(n) for n in range(20, 40)] == TABLE
    for n in range(5, 1001):
        assert evaluate(f, n) == residue_formula(n) == TABLE[n % 20]
        assert residue_formula(n + 20) == residue_formula(n)
    assert [a[n] % 5 for n in range(1, 5)] == [1, 2, 2, 0]
    assert all(c % 5 != 3 for c in a[1:])
    assert 3 not in TABLE
    print('Initial residues n=1..4:', [1, 2, 2, 0])
    print('Candidate period-20 table, indexed by n mod 20:', TABLE)
    print('No residue 3 in independently generated terms n=1..201')
    print('Extraction and period checks pass for n=5..1000')
    print('Scope: exact finite checks; the infinite Hurwitz identification is not a Lean proof')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
