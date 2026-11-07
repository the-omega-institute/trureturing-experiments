#!/usr/bin/env python3
"""Exact, offline checks for Paul D. Hanna's OEIS A091351 conjecture.

Source: https://oeis.org/A091351, published revision 37.
The copied DATA and POWER_ROWS are attributed OEIS data (CC BY-SA 4.0).
All arithmetic below uses Python integers; finite checks are not an unbounded proof.
"""

DATA = [
    1, 1, 1, 1, 2, 1, 1, 4, 3, 1, 1, 9, 9, 4, 1, 1, 24, 30, 16, 5, 1,
    1, 77, 115, 70, 25, 6, 1, 1, 295, 510, 344, 135, 36, 7, 1,
    1, 1329, 2602, 1908, 805, 231, 49, 8, 1,
    1, 6934, 15133, 11904, 5325, 1616, 364, 64, 9, 1,
    1, 41351, 99367, 83028, 39001, 12381, 2919, 540, 81, 10, 1,
]
POWER_ROWS = {
    2: [[1], [2, 1], [4, 4, 1], [9, 14, 6, 1],
        [24, 52, 30, 8, 1], [77, 217, 153, 52, 10, 1]],
    3: [[1], [3, 1], [9, 6, 1], [30, 30, 9, 1],
        [115, 153, 63, 12, 1], [510, 845, 429, 108, 15, 1]],
    4: [[1], [4, 1], [16, 8, 1], [70, 52, 12, 1],
        [344, 336, 108, 16, 1], [1908, 2294, 918, 184, 20, 1]],
}


def zero(n):
    return [[0] * n for _ in range(n)]


def identity(n):
    return [[int(a == b) for b in range(n)] for a in range(n)]


def add(a, b):
    return [[x + y for x, y in zip(r, s)] for r, s in zip(a, b)]


def scale(c, a):
    return [[c * x for x in r] for r in a]


def mul(a, b):
    """Ordinary lower-triangular product; only b <= r <= a can contribute."""
    n = len(a)
    return [[sum(a[i][r] * b[r][j] for r in range(j, i + 1))
             if j <= i else 0 for j in range(n)] for i in range(n)]


def powers(a, k):
    out = [identity(len(a))]
    for _ in range(k):
        out.append(mul(out[-1], a))
    return out


def source_scalar(n):
    """The first FORMULA recurrence, with no quadratic-shift assumption."""
    t = zero(n)
    for a in range(n):
        t[a][0] = 1
        for b in range(1, a + 1):
            t[a][b] = sum(t[a-b][j] * t[j+b-1][b-1]
                          for j in range(a-b+1))
    return t


def quadratic_construction(n):
    u = zero(n)
    for a in range(n):
        u[a][0] = 1
        for b in range(1, a + 1):
            u[a][b] = u[a-1][b-1] + sum(
                u[a-1][r] * u[r-1][b-1] for r in range(b, a))
    return u


def source_rowsums(n):
    """Construct from NAME alone, using only rows strictly before the current row."""
    t = zero(n)
    for a in range(n):
        t[a][0] = 1
        if a:
            earlier = [row[:a] for row in t[:a]]
            pp = powers(earlier, a)
            for b in range(1, a + 1):
                t[a][b] = sum(pp[b][a-b])
    return t


def convolution(a, b):
    return [sum(a[i] * b[k-i] for i in range(k+1)) for k in range(len(a))]


def compose(a, b):
    """Composition modulo z^N; the inner polynomial has constant coefficient zero."""
    assert b[0] == 0
    out = [0] * len(a)
    for c in reversed(a):
        out = convolution(out, b)
        out[0] += c
    return out


def evaluate(coeffs, d):
    out = zero(len(d))
    ident = identity(len(d))
    for c in reversed(coeffs):
        out = add(mul(out, d), scale(c, ident))
    return out


def check_algebra(n, degree_bound=6):
    """Check the draft's intermediate algebra independently, including N=1."""
    u = source_scalar(n)
    ident = identity(n)
    s = [[int(a == b+1) for b in range(n)] for a in range(n)]
    d = mul(s, u)
    assert mul(u, s) == add(d, mul(d, d)), (n, 'quadratic shift')
    assert all(row[0] == 1 for row in u)
    minus_l = add(ident, scale(-1, u))
    inverse_terms = powers(minus_l, n)
    assert inverse_terms[n] == zero(n)
    u_inv = zero(n)
    for term in inverse_terms[:n]:
        u_inv = add(u_inv, term)
    assert mul(u, u_inv) == ident == mul(u_inv, u)
    assert powers(d, n)[n] == zero(n)

    z = [0] * n
    if n > 1:
        z[1] = 1
    g = [0] * n
    if n > 1:
        g[1] = 1
    for r in range(2, n):
        g[r] = -sum(g[i] * g[r-i] for i in range(1, r))
    f = [x + y for x, y in zip(z, convolution(z, z))]
    assert compose(f, g) == z == compose(g, f)
    assert mul(mul(u_inv, d), u) == evaluate(g, d)

    sp = powers(s, degree_bound)
    up = powers(u, degree_bound+1)
    uip = powers(u_inv, degree_bound)
    polynomials = [ident]
    iterate = z
    for m in range(degree_bound+1):
        assert mul(mul(uip[m], d), up[m]) == evaluate(iterate, d)
        assert mul(sp[m], up[m]) == polynomials[m], (n, m, 'factorization')
        # U S^m e = P_m U e, compared as the first columns of the matrices.
        left = mul(u, sp[m])
        right = mul(polynomials[m], u)
        assert [row[0] for row in left] == [row[0] for row in right]
        if m < degree_bound:
            polynomials.append(mul(polynomials[-1], evaluate(iterate, d)))
            iterate = compose(g, iterate)
    for a in polynomials:
        for b in polynomials:
            assert mul(a, b) == mul(b, a)


def main():
    n, max_power = 35, 12
    t = source_scalar(n)
    assert t == quadratic_construction(n)
    # NAME is a third construction, not a restatement of either scalar recurrence.
    assert t == source_rowsums(n)
    flat = [t[a][b] for a in range(11) for b in range(a+1)]
    assert flat == DATA
    print(f'FIDELITY: all {len(DATA)} published DATA terms agree')
    pp = powers(t, n-1)
    fidelity_terms = 0
    for k, rows in POWER_ROWS.items():
        for a, row in enumerate(rows):
            assert pp[k][a][:a+1] == row
            fidelity_terms += len(row)
    print(f'FIDELITY: all {fidelity_terms} displayed entries of powers 2, 3, 4 agree')
    row_checks = 0
    for a in range(n):
        for b in range(a+1):
            assert t[a][b] == sum(pp[b][a-b]), (a, b, 'NAME')
            row_checks += 1
    print(f'DEFINITION: {row_checks} NAME row-sum equalities hold; '
          f'three independent {n} by {n} constructions agree')
    checks = 0
    for p in range(1, max_power+1):
        for q in range(1, max_power+1):
            for r in range(n-max(p, q)+1):
                assert pp[p][r+q-1][q-1] == pp[q][r+p-1][p-1], (r, p, q)
                checks += 1
    print(f'CONJECTURE: {checks} equalities hold for 1 <= p,q <= {max_power} '
          f'and 0 <= n <= {n}-max(p,q)')
    assert pp[2][6][2] == pp[3][5][1] == 845
    print('EXAMPLE: n=4,p=2,q=3 gives (T^2)[6,2]=(T^3)[5,1]=845')
    sizes = (1, 2, 3, 5, 8, 12)
    for size in sizes:
        check_algebra(size)
    print('ROUTE: inverse, quadratic shift, conjugation, factorization, '
          'cyclic-vector and commutation checks pass for N='+','.join(map(str, sizes)))
    print('SCOPE: finite exact arithmetic; no unbounded proof or Lean proof is claimed')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
