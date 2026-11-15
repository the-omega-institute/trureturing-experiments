#!/usr/bin/env python3
"""Exact finite-field checks for Zhu--Gu--Zhang--Li, arXiv:2610.16643v1.

No third-party program is used. Coefficients are ascending, in the prime
fields F_3 and F_5. This checks finite instances, not the universal conjecture.
"""
from itertools import product, combinations
from math import comb

# Literal columns B,C,E,H of sheet1 in DOI 10.5281/zenodo.23279452.
SOURCE_COUNTS = {
    (3, 2, 2): 18, (3, 2, 3): 54, (3, 2, 4): 162, (3, 2, 5): 486,
    (3, 3, 3): 54, (3, 3, 4): 162, (3, 3, 5): 486,
    (3, 3, 6): 1458, (3, 3, 7): 4374, (3, 3, 8): 13122,
    (5, 2, 2): 100, (5, 2, 3): 500, (5, 2, 4): 2500,
    (5, 2, 5): 12500, (5, 2, 6): 62500, (5, 2, 7): 312500,
    (5, 2, 8): 1562500, (5, 2, 9): 7812500,
}


def mul(a, b, q):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] = (c[i + j] + x * y) % q
    return tuple(c)


def split_products(q, cap, degree):
    if degree < 0:
        return []
    out = set()
    for exponents in product(range(cap + 1), repeat=q):
        if sum(exponents) != degree:
            continue
        f = (1,)
        for a, n in enumerate(exponents):
            for _ in range(n):
                f = mul(f, ((-a) % q, 1), q)
        out.add(f)
    return sorted(out)


def dot(a, b, q):
    return sum(x * y for x, y in zip(a, b)) % q


def nullspace(rows, n, q):
    a = [list(row) + [0] * (n - len(row)) for row in rows]
    pivots = []
    r = 0
    for c in range(n):
        pivot = next((i for i in range(r, len(a)) if a[i][c] % q), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        inv = pow(a[r][c] % q, -1, q)
        a[r] = [(x * inv) % q for x in a[r]]
        for i in range(len(a)):
            if i != r:
                z = a[i][c]
                a[i] = [(x - z * y) % q for x, y in zip(a[i], a[r])]
        pivots.append(c)
        r += 1
        if r == len(a):
            break
    basis = []
    for c in range(n):
        if c in pivots:
            continue
        v = [0] * n
        v[c] = 1
        for i, pc in enumerate(pivots):
            v[pc] = -a[i][c] % q
        basis.append(tuple(v))
    assert all(dot(row, v, q) == 0 for row in rows for v in basis)
    return basis


def linear_combinations(basis, n, q):
    for coefficients in product(range(q), repeat=len(basis)):
        yield tuple(sum(c * v[i] for c, v in zip(coefficients, basis)) % q
                    for i in range(n))


def remainder(f, m, q):
    a = list(f)
    for i in range(len(a) - 1, len(m) - 2, -1):
        c = a[i]
        for j, b in enumerate(m):
            a[i - len(m) + 1 + j] = (a[i - len(m) + 1 + j] - c * b) % q
    return a[:len(m) - 1]


def ev(f, q, s):
    return tuple(tuple(sum(comb(n, j) * f[n] * pow(a, n - j, q)
                           for n in range(j, len(f))) % q
                       for a in range(q)) for j in range(s))


def nrt_distance(a, b, q, s):
    return sum(next((s - j for j in range(s) if a[j][x] != b[j][x]), 0)
               for x in range(q))


def fidelity():
    q, s, k = 3, 2, 2
    code = [ev(f, q, s) for f in product(range(q), repeat=k)]
    moduli = split_products(q, s, k + 1)
    distances, deep = [], 0
    images = set()
    for f in product(range(q), repeat=q * s):
        u = ev(f, q, s)
        images.add(u)
        d = min(nrt_distance(u, c, q, s) for c in code)
        distances.append(d)
        criterion = all(remainder(f, m, q)[k] != 0 for m in moduli)
        assert criterion == (d == q * s - k), (f, d, criterion)
        deep += criterion
    assert len(images) == q ** (q * s)  # Hermite interpolation bijectivity.
    assert max(distances) == 4 and deep == 18
    print('FIDELITY q=3 s=2 k=2: all 729 words; radius=4; deep=18; '
          'Hermite bijection and Theorems 3.1/3.3 agree with direct NRT distances')


def small_block(q, degree):
    """Enumerate every functional on P_{<=degree}; no Segre input is used."""
    qs = split_products(q, 1, degree)
    return {l for l in product(range(q), repeat=degree + 1)
            if all(dot(l, p, q) for p in qs)}


def classification(q, s, k, quadratic_block):
    t = q * s - k - 1
    qs = split_products(q, s, t)
    forced = []
    # For R capped at s-1, every (X-a)R is admissible. Exhaustively
    # enumerating the q affine tests excludes L(R) != 0, forcing L(R)=0.
    for r in split_products(q, s - 1, t - 1):
        for a in range(q):
            assert mul(r, ((-a) % q, 1), q) in qs
        forced.append(r)
    # For q=5, enumerate ALL 125 functionals on degree <=2 once.
    # Every restriction P -> L(RP) must be in that computed set. The
    # computed set forces L(R)=L(XR)=0. These are necessary consequences
    # of the original constraints, not assumptions of the conjecture.
    if q == 5 and t >= 2:
        assert quadratic_block == {(0, 0, c) for c in range(1, q)}
        squarefree = split_products(q, 1, 2)
        for r in split_products(q, s - 1, t - 2):
            assert all(mul(r, p, q) in qs for p in squarefree)
            forced.extend([r, (0,) + r])
    basis = nullspace(forced, t + 1, q)
    solutions = set()
    count = 0
    for l in linear_combinations(basis, t + 1, q):
        count += 1
        if all(dot(l, p, q) for p in qs):
            solutions.add(l)
    expected = {(0,) * t + (c,) for c in range(1, q)}
    assert solutions == expected, (q, s, k, solutions, expected)
    # Complement coefficients give a triangular isomorphism from the
    # high coefficients f_k,...,f_{N-1} to these functionals. Check the
    # actual source remainder matrix, its invertibility, and that the
    # surviving vectors have exactly f_k != 0 and all higher f_i = 0.
    n = q * s
    moduli = split_products(q, s, k + 1)
    pi_s = split_products(q, s, n)
    assert len(pi_s) == 1
    matrix = []
    for m in moduli:
        row = tuple(remainder((0,) * i + (1,), m, q)[k] for i in range(k, n))
        quotient = [0] * (t + 1)
        a = list(pi_s[0])
        for i in range(n, k, -1):
            c = a[i]
            quotient[i - (k + 1)] = c
            for j, b in enumerate(m):
                a[i - (k + 1) + j] = (a[i - (k + 1) + j] - c * b) % q
        assert not any(a)
        matrix.append((tuple(quotient), row))
    # Solve for the universal pairing matrix T with row = Q*T.
    # Recover the formal-series pairing directly, then check every
    # source modulus. Split products need not span at very high t.
    denominator = pi_s[0]
    inv = [0] * n
    inv[0] = 1
    for j in range(1, n):
        inv[j] = -sum(denominator[n - h] * inv[j - h]
                      for h in range(1, min(j, n) + 1)) % q
    transform = [tuple(inv[j + i + k - n + 1] if 0 <= j + i + k - n + 1 < n else 0
                       for i in range(t + 1)) for j in range(t + 1)]
    assert len(nullspace(transform, t + 1, q)) == 0
    for a, b in matrix:
        assert tuple(sum(a[j] * transform[j][i] for j in range(t + 1)) % q
                     for i in range(t + 1)) == b
    # T*[c,0,...,0] is exactly [0,...,0,c].
    assert tuple(row[0] for row in transform) == (0,) * t + (1,)
    deep_words = len(solutions) * q ** k
    print(f'CLASSIFICATION q={q} s={s} k={k}: profiles={len(qs)} '
          f'residual_dimension={len(basis)} enumerated={count} '
          f'deep_cosets={len(solutions)} deep_words={deep_words}')
    return deep_words


def main():
    fidelity()
    block = small_block(5, 2)
    assert block == {(0, 0, c) for c in range(1, 5)}
    print('SQUAREFREE_BLOCK q=5 degree=2: all 125 functionals, 4 survivors')
    comparisons = 0
    for q, s in [(3, 2), (3, 3), (5, 2), (5, 3)]:
        for k in range(max(2, s), q * s):
            count = classification(q, s, k, block)
            # DOI 10.5281/zenodo.23279452, NRT deep hole.xlsx, sheet1:
            # all 18 corresponding rows for (3,2), (3,3), (5,2) have
            # H=(q-1)*q^k. There are no (5,3) rows in that workbook.
            if (q, s, k) in SOURCE_COUNTS:
                assert count == SOURCE_COUNTS[q, s, k]
                comparisons += 1
    assert comparisons == 18
    print('SOURCE_DATA: 18/18 overlapping workbook rows agree; '
          '12 q=5,s=3 cases extend the workbook; 30 cases total')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
