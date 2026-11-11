#!/usr/bin/env python3
"""Exact arithmetic diagnostics for Orzel's higher-characteristic converse.

No external dependencies. F9 = F3[i]/(i^2+1); a+bi is encoded by a+3*b.
The geometric nonmembership argument is in the preregistration draft. Finite
point checks below do not establish nonmembership over an algebraic closure.
"""

from collections import Counter
from itertools import product
import json
import random


def fadd(x, y):
    return (x % 3 + y % 3) % 3 + 3 * ((x // 3 + y // 3) % 3)


def fmul(x, y):
    a, b, c, d = x % 3, x // 3, y % 3, y // 3
    return (a * c + 2 * b * d) % 3 + 3 * ((a * d + b * c) % 3)


ADD = [[fadd(x, y) for y in range(9)] for x in range(9)]
MUL = [[fmul(x, y) for y in range(9)] for x in range(9)]
NEG = [next(y for y in range(9) if ADD[x][y] == 0) for x in range(9)]
INV = [0] + [next(y for y in range(1, 9) if MUL[x][y] == 1) for x in range(1, 9)]


def fpow(x, k):
    out = 1
    for _ in range(k):
        out = MUL[out][x]
    return out


def const(n, c):
    return {(0,) * n: c} if c else {}


def var(n, i):
    return {tuple(int(j == i) for j in range(n)): 1}


def add(*polys):
    out = {}
    for p in polys:
        for m, c in p.items():
            out[m] = ADD[out.get(m, 0)][c]
    return {m: c for m, c in out.items() if c}


def scale(c, p):
    return {m: MUL[c][v] for m, v in p.items() if MUL[c][v]}


def neg(p):
    return {m: NEG[c] for m, c in p.items()}


def sub(p, q):
    return add(p, neg(q))


def mul(p, q):
    out = {}
    for m, c in p.items():
        for t, v in q.items():
            u = tuple(a + b for a, b in zip(m, t, strict=True))
            out[u] = ADD[out.get(u, 0)][MUL[c][v]]
    return {m: c for m, c in out.items() if c}


def power(p, k, n):
    out = const(n, 1)
    for _ in range(k):
        out = mul(out, p)
    return out


def esym(forms, degree, n):
    out = [const(n, 1)] + [{} for _ in range(degree)]
    for form in forms:
        for d in range(degree, 0, -1):
            out[d] = add(out[d], mul(form, out[d - 1]))
    return out


def deriv(p, j):
    out = {}
    for m, c in p.items():
        if m[j] % 3:
            t = list(m)
            t[j] -= 1
            out[tuple(t)] = MUL[c][m[j] % 3]
    return out


def hessian(p, n):
    return [[deriv(deriv(p, i), j) for j in range(n)] for i in range(n)]


def evaluate(p, point):
    out = 0
    for m, c in p.items():
        term = c
        for x, k in zip(point, m, strict=True):
            term = MUL[term][fpow(x, k)]
        out = ADD[out][term]
    return out


def rank3(matrix):
    a = [row[:] for row in matrix]
    r = 0
    for j in range(len(a[0])):
        pivot = next((i for i in range(r, len(a)) if a[i][j] % 3), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        a[r] = [(x * a[r][j]) % 3 for x in a[r]]
        for i in range(r + 1, len(a)):
            c = a[i][j]
            if c:
                a[i] = [(x - c * y) % 3 for x, y in zip(a[i], a[r], strict=True)]
        r += 1
        if r == len(a):
            break
    return r


def rank9(matrix):
    a = [row[:] for row in matrix]
    r = 0
    for j in range(len(a[0])):
        pivot = next((i for i in range(r, len(a)) if a[i][j]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        c = INV[a[r][j]]
        a[r] = [MUL[x][c] for x in a[r]]
        for i in range(r + 1, len(a)):
            c = a[i][j]
            if c:
                a[i] = [ADD[x][NEG[MUL[c][y]]] for x, y in zip(a[i], a[r], strict=True)]
        r += 1
        if r == len(a):
            break
    return r


def linear(n, rng):
    return add(*(scale(rng.randrange(9), var(n, j)) for j in range(n)))


def family(a, b):
    n = a + b
    x = [var(n, j) for j in range(a)]
    y = [var(n, a + j) for j in range(b)]
    aa = add(*(mul(z, z) for z in x))
    bb = add(*(mul(z, z) for z in y))
    return x, y, aa, bb, mul(aa, bb)


def powers_block(form):
    # alpha=1+i, alpha^4=-1; all roots of z^4+1 are alpha*i^j.
    roots = [MUL[4][fpow(3, j)] for j in range(4)]
    assert len(set(roots)) == 4 and all(fpow(x, 4) == 2 for x in roots)
    return [scale(NEG[x], form) for x in roots]


def square_representation(diagonal_forms, n):
    # e4(L1,-L1,...) = sum_{i<j} Li^2 Lj^2; alpha^4=-1 negates it.
    forms = [scale(c, form) for form in diagonal_forms for c in (4, NEG[4])]
    for form in diagonal_forms:
        forms.extend(powers_block(form))
    q = add(*(mul(form, form) for form in diagonal_forms))
    assert esym(forms, 4, n)[4] == mul(q, q)
    return forms, mul(q, q)


def source_fidelity():
    rng = random.Random(250905009)
    n = 4
    for m in range(11):
        forms = [linear(n, rng) for _ in range(m)]
        e = esym(forms, 4, n)
        ps = [{}] + [add(*(power(form, d, n) for form in forms)) for d in range(1, 5)]
        # Appendix A, Newton identities: 4 e4 = e3 p1 - e2 p2 + e1 p3 - p4.
        rhs = add(mul(e[3], ps[1]), neg(mul(e[2], ps[2])), mul(e[1], ps[3]), neg(ps[4]))
        assert e[4] == rhs
        normal = add(mul(e[1], sub(e[3], mul(e[1], e[2]))),
                     neg(mul(e[2], e[2])), power(e[1], 4, n), neg(ps[4]))
        assert e[4] == normal
        # Lemma 3.1: freely append one fourth power to a Sym form.
        extra = linear(n, rng)
        block = powers_block(extra)
        be = esym(block, 4, n)
        assert be[1:4] == [{}, {}, {}] and be[4] == power(extra, 4, n)
        assert esym(forms + block, 4, n)[4] == add(e[4], power(extra, 4, n))
    print('FIDELITY Newton degree-4 identity and Lemma 3.1: m=0..10 over F9 OK')


def generic_hessian_identity():
    rng = random.Random(307)
    n = 4
    ell = linear(n, rng)
    q = add(*(mul(linear(n, rng), linear(n, rng)) for _ in range(3)))
    c = mul(mul(linear(n, rng), linear(n, rng)), linear(n, rng))
    w = add(*(power(linear(n, rng), 4, n) for _ in range(3)))
    f = add(mul(ell, c), neg(mul(q, q)), w)
    hf, hc, hq, hw = (hessian(p, n) for p in (f, c, q, w))
    for i in range(n):
        for j in range(n):
            assert hw[i][j] == {}
            expected = add(mul(ell, hc[i][j]), scale(NEG[2], mul(q, hq[i][j])),
                           mul(deriv(ell, i), deriv(c, j)), mul(deriv(c, i), deriv(ell, j)),
                           scale(NEG[2], mul(deriv(q, i), deriv(q, j))))
            assert hf[i][j] == expected
    print('NORMAL_FORM exact symbolic Hessian product rule and H(fourth powers)=0 OK')


def symbolic_family_checks():
    for a, b in ((4, 4), (4, 5), (5, 4), (5, 5), (4, 6), (6, 4)):
        n = a + b
        x, y, aa, bb, f = family(a, b)
        h = hessian(f, n)
        for i in range(n):
            for j in range(n):
                if i < a and j < a:
                    expected = scale(2, bb) if i == j else {}
                elif i >= a and j >= a:
                    expected = scale(2, aa) if i == j else {}
                else:
                    expected = mul(x[i] if i < a else x[j], y[j - a] if j >= a else y[i - a])
                assert h[i][j] == expected
        first = mul(aa, add(mul(y[0], y[0]), mul(y[1], y[1])))
        second = mul(aa, add(*(mul(z, z) for z in y[2:])))
        assert first and second and add(first, second) == f
        assert sub(mul(add(aa, bb), add(aa, bb)), mul(sub(aa, bb), sub(aa, bb))) == f
    print('FAMILY symbolic Hessians, two nonzero reducible summands, square difference: six (a,b) pairs OK')


def finite_rank_checks():
    for a, b in ((4, 4), (4, 5), (5, 4)):
        n = a + b
        counts = Counter()
        for point in product(range(3), repeat=n):
            x, y = point[:a], point[a:]
            aa, bb = sum(v * v for v in x) % 3, sum(v * v for v in y) % 3
            h = [[(2 * bb * int(i == j)) % 3 if i < a and j < a else
                  (2 * aa * int(i == j)) % 3 if i >= a and j >= a else
                  (x[i] * y[j - a]) % 3 if i < a else (y[i - a] * x[j]) % 3
                  for j in range(n)] for i in range(n)]
            r = rank3(h)
            assert (r <= 3) == (aa == 0 and bb == 0)
            if aa and bb:
                assert r == n - 1
            elif aa == 0 and bb:
                assert r == a
            elif aa and bb == 0:
                assert r == b
            else:
                assert r in (0, 2)
            counts[str((bool(aa), bool(bb), r))] += 1
        print(f'RANK_F3 a={a} b={b} points={3**n} strata=' + json.dumps(dict(sorted(counts.items()))))
    rng = random.Random(904)
    for a, b in ((4, 4), (4, 5), (5, 4)):
        n = a + b
        x, y, aa, bb, f = family(a, b)
        h = hessian(f, n)
        for _ in range(200):
            point = [rng.randrange(9) for _ in range(n)]
            av, bv = evaluate(aa, point), evaluate(bb, point)
            assert (rank9([[evaluate(p, point) for p in row] for row in h]) <= 3) == (av == bv == 0)
    print('RANK_F9 600 seeded extension-field points OK (finite diagnostics only)')


def nullcone_and_upper_bound():
    for d in range(4, 7):
        vectors = []
        for j in range(1, d):
            for c in (3, NEG[3]):
                v = [0] * d
                v[0], v[j] = 1, c
                assert ADD[fmul(v[0], v[0])][fmul(v[j], v[j])] == 0
                vectors.append(v)
        assert rank9(vectors) == d
    print('NULLCONES explicit isotropic spanning vectors in dimensions 4..6 over F9 OK')
    for a, b in ((4, 4), (4, 5), (5, 4)):
        n = a + b
        x, y, aa, bb, f = family(a, b)
        plus_forms, plus = square_representation(x + y, n)
        minus_forms, minus = square_representation(x + [scale(3, z) for z in y], n)
        assert sub(plus, minus) == f
        print(f'SYM_UPPER_BOUND a={a} b={b}: AB=e4(list_{len(plus_forms)})-e4(list_{len(minus_forms)}) OK')


def main():
    assert fpow(3, 2) == 2 and fpow(4, 4) == 2
    assert all(MUL[x][INV[x]] == 1 for x in range(1, 9))
    source_fidelity()
    generic_hessian_identity()
    symbolic_family_checks()
    finite_rank_checks()
    nullcone_and_upper_bound()
    print('SCOPE: nonmembership and minimum top fan-in require the geometric argument; finite checks alone do not prove them.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
