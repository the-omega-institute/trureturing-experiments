#!/usr/bin/env python3
"""Exact checks for Wu, arXiv:2607.23182v1, Appendix L.7.

Python standard library only. All arithmetic is integral or Fraction arithmetic.
The finite tests do not replace the proposed argument for arbitrary dimension.
"""
from fractions import Fraction as Q
from itertools import product


def vector(xs):
    return tuple(map(Q, xs))


def matrix(rows):
    return tuple(vector(row) for row in rows)


def eye(d):
    return matrix([[int(i == j) for j in range(d)] for i in range(d)])


def add(x, y):
    return tuple(a + b for a, b in zip(x, y))


def sub(x, y):
    return tuple(a - b for a, b in zip(x, y))


def scale(a, x):
    return tuple(a * b for b in x)


def dot(x, y):
    return sum((a * b for a, b in zip(x, y)), Q(0))


def mv(a, x):
    return tuple(dot(row, x) for row in a)


def mm(a, b):
    return tuple(tuple(dot(row, col) for col in zip(*b)) for row in a)


def columns(a):
    return tuple(zip(*a))


def from_columns(cols):
    return tuple(tuple(row) for row in zip(*cols))


def determinant(a):
    if not a:
        return Q(1)
    return sum(((-1) ** j * a[0][j] * determinant(tuple(
        tuple(row[k] for k in range(len(a)) if k != j) for row in a[1:]))
        for j in range(len(a))), Q(0))


def rank(a):
    if not a:
        return 0
    a = [list(row) for row in a]
    r = 0
    for j in range(len(a[0])):
        pivot = next((i for i in range(r, len(a)) if a[i][j]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        a[r] = [x / a[r][j] for x in a[r]]
        for i in range(len(a)):
            if i != r:
                c = a[i][j]
                a[i] = [x - c * y for x, y in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def inverse(a):
    d = len(a)
    rows = [list(row) + list(unit) for row, unit in zip(a, eye(d))]
    for j in range(d):
        pivot = next(i for i in range(j, d) if rows[i][j])
        rows[j], rows[pivot] = rows[pivot], rows[j]
        rows[j] = [x / rows[j][j] for x in rows[j]]
        for i in range(d):
            if i != j:
                c = rows[i][j]
                rows[i] = [x - c * y for x, y in zip(rows[i], rows[j])]
    inv = tuple(tuple(row[d:]) for row in rows)
    assert mm(a, inv) == eye(d) == mm(inv, a)
    return inv


def normal(b):
    """Coefficients of ell(v) = det[v,B] by expansion in the first column."""
    d = len(b)
    return tuple((-1) ** i * determinant(tuple(row for j, row in enumerate(b) if j != i))
                 for i in range(d))


def in_slab(i, t, interior=False):
    if interior:
        return t < 0 if i == 0 else (0 < t < 1 if i == 1 else 1 < t)
    return t <= 0 if i == 0 else (0 <= t <= 1 if i == 1 else 1 <= t)


def continuity_dimension(d):
    """Rank of ALL continuity equations on three arbitrary affine rules."""
    width = 3 * (d * d + d)
    def pos(i, r, c):
        return i * (d * d + d) + r * (d + 1) + c
    equations = []
    for r in range(d):
        for c in range(1, d + 1):
            row = [Q(0)] * width
            row[pos(0, r, c)], row[pos(1, r, c)] = Q(1), Q(-1)
            equations.append(row)
        for c in range(1, d):
            row = [Q(0)] * width
            row[pos(1, r, c)], row[pos(2, r, c)] = Q(1), Q(-1)
            equations.append(row)
        row = [Q(0)] * width
        row[pos(1, r, 0)] = row[pos(1, r, d)] = Q(1)
        row[pos(2, r, 0)] = row[pos(2, r, d)] = Q(-1)
        equations.append(row)
    assert rank(equations) == 2 * d * d
    assert width - rank(equations) == d * d + 3 * d
    return width - rank(equations)


def check_family(a, offsets, label):
    d = len(a[0])
    b = tuple(tuple(row[j] for j in range(1, d)) for row in a[0])
    assert all(tuple(tuple(row[j] for j in range(1, d)) for row in ai) == b for ai in a)
    first = [columns(ai)[0] for ai in a]
    assert offsets[0] == offsets[1]
    assert add(first[1], offsets[1]) == add(first[2], offsets[2])
    # These coefficient equalities certify continuity for EVERY transverse z.
    ell = normal(b)
    s = tuple(determinant(ai) for ai in a)
    assert s == tuple(dot(ell, col) for col in first)
    assert s[0] > 0 > s[1] and s[2] > 0
    assert rank(b) == d - 1 and any(ell)
    assert all(dot(ell, col) == 0 for col in columns(b))
    assert len(set(zip(a, offsets))) == 3
    inverses = tuple(inverse(ai) for ai in a)
    # Identity of affine functionals for ALL y, not a sampled inverse formula.
    for i, inv in enumerate(inverses):
        assert inv[0] == scale(Q(1) / s[i], ell)
        shift = -dot(inv[0], sub(offsets[i], offsets[0]))
        assert shift == (Q(0) if i < 2 else Q(1) - s[1] / s[2])
    # Thus inverse t coordinates are u/s1, u/s2, 1+(u-s2)/s3.
    # Their slab inequalities give u<=0, s2<=u<=0, u>=s2.
    us = (s[1] - 1, s[1], s[1] / 2, Q(0), Q(1))
    transverse_values = (Q(-2, 3), Q(0), Q(5, 2))
    target_count = 0
    for u in us:
        for z in product(transverse_values, repeat=d - 1):
            y = add(offsets[0], add(scale(u / s[0], first[0]), mv(b, z)))
            assert dot(ell, sub(y, offsets[0])) == u
            xs = [mv(inv, sub(y, bi)) for inv, bi in zip(inverses, offsets)]
            active = {i for i, x in enumerate(xs) if in_slab(i, x[0])}
            expected = ({0} if u <= 0 else set()) | ({1} if s[1] <= u <= 0 else set())
            expected |= {2} if u >= s[1] else set()
            assert active == expected
            interior_active = {i for i, x in enumerate(xs) if in_slab(i, x[0], True)}
            if u not in (s[1], 0):
                assert interior_active == ({0} if u < s[1] else {0, 1, 2} if u < 0 else {2})
            preimages = {xs[i] for i in active}
            assert len(preimages) == (1 if u < s[1] or u > 0 else 3 if s[1] < u < 0 else 2)
            for i in active:
                assert add(mv(a[i], xs[i]), offsets[i]) == y
            target_count += 1
    chambers = ({0}, {0, 1, 2}, {2})
    differences = tuple((len(c - e), len(e - c)) for c, e in zip(chambers, chambers[1:]))
    assert differences == ((0, 2), (2, 0))
    assert all(1 not in pair for pair in differences)
    print(f'{label}: d={d}; determinants={s}; exact targets={target_count}; '
          f'chamber fibers=1,3,1; adjacent directional differences={differences}')
    return target_count


def generated(d, k):
    t = [list(row) for row in eye(d)]
    # Rational elementary row operations, each preserving determinant 1.
    for i in range(d):
        j = (i + 1) % d
        t[i] = [x + Q(k + i + 1, k + 2) * y for x, y in zip(t[i], t[j])]
    t = matrix(t)
    assert determinant(t) == 1
    cols = columns(t)
    b = from_columns(cols[1:])
    slopes = (Q(k + 1, k + 2), -Q(k + 2, k + 3), Q(k + 3, k + 4))
    first = [add(scale(slopes[i], cols[0]),
                 mv(b, vector(Q((i + 1) * (j + 1), k + 2) for j in range(d - 1))))
             for i in range(3)]
    a = tuple(from_columns((ai,) + cols[1:]) for ai in first)
    offset = vector(Q(j - k, k + 3) for j in range(d))
    return a, (offset, offset, add(offset, sub(first[1], first[2])))


def main():
    # Fidelity: exactly Example A of Appendix L.3, including its breakpoint values.
    a = (matrix([[1]]), matrix([[-1]]), matrix([[1]]))
    offsets = (vector([0]), vector([0]), vector([-2]))
    check_family(a, offsets, 'FIDELITY Example A')
    def source_f(x):
        return x if x <= 0 else -x if x <= 1 else x - 2
    assert source_f(Q(0)) == 0 and source_f(Q(1)) == -1
    for y, expected in [(Q(-2), (Q(-2),)), (Q(-1, 2), (Q(-1, 2), Q(1, 2), Q(3, 2))),
                        (Q(1), (Q(3),))]:
        candidates = (y, -y, y + 2)
        valid = tuple(x for i, x in enumerate(candidates)
                      if (x <= 0 if i == 0 else 0 < x <= 1 if i == 1 else x > 1))
        assert valid == expected and all(source_f(x) == y for x in valid)
    print('FIDELITY: Example A formula, continuity, branch images, fibers and SB failure checked')
    a2 = (matrix([[2, 1], [1, 1]]), matrix([[1, 1], [2, 1]]),
          matrix([[3, 1], [2, 1]]))
    b2 = (vector([0, 0]), vector([0, 0]), vector([-2, 0]))
    total = check_family(a2, b2, 'EXPLICIT rational pair')
    # Explicit fully mixed 3D example, not merely the coordinate product witness.
    a3 = (matrix([[2, 1, 0], [1, 1, 1], [1, 0, 2]]),
          matrix([[-1, 1, 0], [1, 1, 1], [0, 0, 2]]),
          matrix([[4, 1, 0], [1, 1, 1], [6, 0, 2]]))
    b3 = (vector([Q(1, 2), -1, 2]), vector([Q(1, 2), -1, 2]),
          vector([Q(-9, 2), -1, -4]))
    total += check_family(a3, b3, 'EXPLICIT rational triple')
    for d in (2, 3):
        for k in range(8):
            total += check_family(*generated(d, k), f'FAMILY k={k}')
        dim = continuity_dimension(d)
        print(f'CONTINUITY d={d}: constraint rank={2*d*d}; full compatible dimension={dim}')
    print(f'SETTLEMENT CHECK: 18 exact rational maps in dimensions 2 and 3; {total} targets')
    print('GENERAL SCOPE: finite checks and affine identities; arbitrary-d existence awaits Lean proof')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
