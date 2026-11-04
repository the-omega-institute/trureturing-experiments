#!/usr/bin/env python3
"""Exact finite checks for Clarke, arXiv:2610.17988v1, Conjecture 8.7.

Only integer/rational arithmetic contributes evidence. These are finite checks,
not a proof of the conjecture or of the absence of later literature.
"""
from itertools import combinations, product
from math import gcd, lcm
from random import Random
import hashlib
import json
from pathlib import Path
import sympy as sp

ZERO, ONE, EPS = (0, 0), (1, 0), (0, 1)
SEED, MAX_M, MAX_DEGREE = 261087, 7, 8


def add(a, b):
    return a[0] + b[0], a[1] + b[1]


def mul(a, b):
    return a[0] * b[0] + a[1] * b[1], a[0] * b[1] + a[1] * b[0]


def convolution(a, b, bound=None):
    size = len(a) + len(b) - 1
    if bound is not None:
        size = min(size, bound + 1)
    out = [ZERO] * size
    for i, x in enumerate(a):
        for j, y in enumerate(b[:max(0, size - i)]):
            out[i + j] = add(out[i + j], mul(x, y))
    return out


def power(a, n):
    out = [ONE]
    for _ in range(n):
        out = convolution(out, a)
    return out


def denominator(p, q):
    return convolution(power([ONE, (-1, 0)], p), power([ONE, (0, -1)], q))


def kernel(p, q, m):
    assert m > 0 and m % 2 == 1
    return convolution(power([ONE] * m, p),
                       power([ONE if j % 2 == 0 else EPS for j in range(m)], q))


def coeff(a, k):
    return a[k] if k < len(a) else ZERO


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def affine_action(A, b, x, k=1):
    return tuple(dot(row, x) + k * z for row, z in zip(A, b))


def facets(vertices):
    """Enumerate supporting hyperplanes with exact nullspaces; dimension <= 3."""
    n = len(vertices[0])
    assert sp.Matrix([[x - y for x, y in zip(v, vertices[0])]
                      for v in vertices[1:]]).rank() == n
    result = set()
    for face in combinations(vertices, n):
        if n == 1:
            normal = [1]
        else:
            matrix = sp.Matrix([[x - y for x, y in zip(v, face[0])] for v in face[1:]])
            null = matrix.nullspace()
            if len(null) != 1:
                continue
            scale = lcm(*(int(x.q) for x in null[0]))
            normal = [int(x * scale) for x in null[0]]
        c = dot(normal, face[0])
        differences = [dot(normal, v) - c for v in vertices]
        if all(d <= 0 for d in differences):
            pass
        elif all(d >= 0 for d in differences):
            normal, c = [-x for x in normal], -c
        else:
            continue
        divisor = gcd(*(abs(x) for x in normal))
        result.add((tuple(x // divisor for x in normal), c // divisor))
    assert result
    return sorted(result)


def exact_counts(vertices, inequalities, A, b, k):
    """Count integer fibers of kP using exact inclusive integer intervals.

    The last coordinate is eliminated using all facet inequalities; fixed-point
    equations either determine that coordinate or leave the whole interval.
    """
    n = len(vertices[0])
    if k == 0:
        return 1, 1
    lower = [k * min(v[i] for v in vertices) for i in range(n)]
    upper = [k * max(v[i] for v in vertices) for i in range(n)]
    total = fixed = 0
    for prefix in product(*(range(lower[i], upper[i] + 1) for i in range(n - 1))):
        lo, hi = lower[-1], upper[-1]
        for normal, c in inequalities:
            last = normal[-1]
            rhs = k * c - dot(normal[:-1], prefix)
            if last > 0:
                hi = min(hi, rhs // last)
            elif last < 0:
                lo = max(lo, -((-rhs) // last))
            elif rhs < 0:
                hi = lo - 1
                break
        if lo > hi:
            continue
        total += hi - lo + 1
        unique, possible = None, True
        for i in range(n):
            coefficient = A[i][-1] - int(i == n - 1)
            constant = dot(A[i][:-1], prefix) + k * b[i]
            if i < n - 1:
                constant -= prefix[i]
            if coefficient == 0:
                if constant:
                    possible = False
                    break
            elif (-constant) % coefficient:
                possible = False
                break
            else:
                value = (-constant) // coefficient
                if unique is not None and unique != value:
                    possible = False
                    break
                unique = value
        if possible:
            fixed += hi - lo + 1 if unique is None else int(lo <= unique <= hi)
    assert 0 <= fixed <= total and (total - fixed) % 2 == 0
    return total, fixed


def brute_counts(vertices, inequalities, A, b, k):
    n = len(vertices[0])
    if k == 0:
        return 1, 1
    points = [x for x in product(*(range(k * min(v[i] for v in vertices),
                                      k * max(v[i] for v in vertices) + 1)
                                  for i in range(n)))
              if all(dot(a, x) <= k * c for a, c in inequalities)]
    return len(points), sum(affine_action(A, b, x, k) == x for x in points)


def character(N, F):
    assert (N + F) % 2 == 0 and (N - F) % 2 == 0
    return (N + F) // 2, (N - F) // 2


def source_examples():
    square = list(product(range(2), repeat=2))
    pyramid = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 0, 1), (0, 1, 1)]
    bipyramid = [(x, y, 0) for x, y in square] + [(0, 0, 1), (1, 1, -1)]
    tetrahedron = [(0, 0, 0), (2, 0, 0), (1, 1, 0), (1, 0, 1)]
    txi = [v + (z,) for v in tetrahedron for z in range(2)]
    return [
        ('source-Example-2.2-orthogonal-square', square, ((-1, 0), (0, 1)), (1, 0)),
        ('source-Example-2.3-diagonal-square', square, ((0, 1), (1, 0)), (0, 0)),
        ('source-Example-5.1-pyramid', pyramid, ((0, 1, 0), (1, 0, 0), (0, 0, 1)), (0, 0, 0)),
        ('source-Example-5.2-bipyramid', bipyramid, ((-1, 0, -1), (0, 1, 0), (0, 0, 1)), (1, 0, 0)),
        ('source-Example-5.3-TxI', txi,
         ((1, 0, 0, 0), (0, 0, 1, 0), (0, 1, 0, 0), (0, 0, 0, 1)), (0, 0, 0, 0)),
    ]


def txi_counts(k):
    """Enumerate the exact integer fibers of T x [0,1].

    T has y1,y2 >= 0, y1+y2 <= k, and y1+y2 <= x <= 2k-y1-y2.
    Fixed points have y1=y2; the fourth coordinate ranges from 0 to k.
    """
    total = (k + 1) * sum((s + 1) * (2 * (k - s) + 1) for s in range(k + 1))
    fixed = (k + 1) * sum(2 * k - 4 * a + 1 for a in range(k // 2 + 1))
    return total, fixed


def random_examples():
    rng = Random(SEED)
    out = []
    for n in range(1, 4):
        for mode in ['central', 'reflection', 'affine-reflection']:
            for sample in range(3):
                A = tuple(tuple((-1 if mode == 'central' or i == 0 else 1)
                                if i == j else 0 for j in range(n)) for i in range(n))
                b = (1,) + (0,) * (n - 1) if mode == 'affine-reflection' else (0,) * n
                for _ in range(100):
                    points = [tuple(rng.randrange(-1, 2) for _ in range(n))
                              for _ in range(n + 3)]
                    vertices = sorted(set(points + [affine_action(A, b, x) for x in points]))
                    if sp.Matrix([[x - y for x, y in zip(v, vertices[0])]
                                  for v in vertices[1:]]).rank() == n:
                        break
                else:
                    raise AssertionError('Failed to generate a full-dimensional polytope')
                out.append((f'random-d{n}-{mode}-{sample}', vertices, A, b))
    return out


def main():
    print('Clarke arXiv:2610.17988v1 Conjecture 8.7: exact finite checks')
    print('Python integer arithmetic; SymPy', sp.__version__, 'exact rank/nullspace only')
    print('seed', SEED, 'dimensions 1..3; dilations 1..7; numerator degrees 0..8')
    count = 0
    for p, q, m in product(range(1, 6), range(1, 6), range(1, 10, 2)):
        K = kernel(p, q, m)
        assert all(min(x) >= 0 for x in K)
        lhs = convolution(denominator(p, q), K)
        D = denominator(p, q)
        rhs = [ZERO] * (m * (len(D) - 1) + 1)
        for j, a in enumerate(D):
            rhs[m * j] = a
        assert lhs == rhs
        count += 1
    print('D(t) K_m(t) = D(t^m), effective K_m:', count, 'exact polynomial identities')
    all_readings = []
    effective_premises = controls = odd_identities = brute_crosschecks = 0
    for name, vertices, A, b in source_examples() + random_examples():
        n = len(vertices[0])
        assert set(affine_action(A, b, x) for x in vertices) == set(vertices)
        assert all(affine_action(A, b, affine_action(A, b, x)) == x for x in vertices)
        inequalities = facets(vertices)
        homogeneous = sp.Matrix([list(row) + [z] for row, z in zip(A, b)]
                                + [[0] * n + [1]])
        p = len((homogeneous - sp.eye(n + 1)).nullspace())
        q = len((homogeneous + sp.eye(n + 1)).nullspace())
        assert p + q == n + 1 and p >= 1 and q >= 1
        counts = [txi_counts(k) if name.endswith('TxI') else
                  exact_counts(vertices, inequalities, A, b, k)
                  for k in range(MAX_M * MAX_DEGREE + 1)]
        for k in range(4):
            assert counts[k] == brute_counts(vertices, inequalities, A, b, k)
            brute_crosschecks += 1
        E = [character(N, F) for N, F in counts]
        D = denominator(p, q)
        phi = convolution(D, E, MAX_M * MAX_DEGREE)
        if name.endswith('orthogonal-square'):
            assert [a - b for a, b in phi[:9]] == [1, -1, 2, -2, 2, -2, 2, -2, 2]
            assert [a + b for a, b in phi[:9]] == [1, 1, 0, 0, 0, 0, 0, 0, 0]
            print('FIDELITY Example 2.2: phi(sigma) = (1+t^2)/(1+t), degrees 0..8 reproduced')
        if name.endswith('diagonal-square'):
            assert phi == [ONE, ONE] + [ZERO] * (len(phi) - 2)
            print('FIDELITY Example 2.3: phi = (1+t)*trivial, degrees 0..56 reproduced')
        if name.endswith('5.1-pyramid'):
            assert all(F == (k // 2 + 1) ** 2 for k, (_, F) in enumerate(counts))
            assert [a + b for a, b in phi] == [1, 1] + [0] * (len(phi) - 2)
            assert [a - b for a, b in phi[:9]] == [1, -1, 2, -2, 2, -2, 2, -2, 2]
            print('FIDELITY Example 5.1: h* = 1+t; N(k) = (floor(k/2)+1)^2, k=0..56')
        if name.endswith('5.2-bipyramid'):
            assert phi == [ONE, (1, 1), ONE] + [ZERO] * (len(phi) - 3)
            print('FIDELITY Example 5.2: phi = 1+(trivial+sign)*t+t^2, degrees 0..56')
        if name.endswith('TxI'):
            assert all(2 * F == (k + 1) ** 2 * (k + 2) for k, (_, F) in enumerate(counts))
            assert phi == [ONE, (4, 1), (2, 0)] + [ZERO] * (len(phi) - 3)
            print('FIDELITY Example 5.3: phi=(1+4t+2t^2)*trivial+t*sign; obstructing-face example')
        premise = all(min(x) >= 0 for x in phi)
        effective_premises += premise
        controls += not premise
        outcomes = []
        for m in range(1, MAX_M + 1):
            actual = convolution(D, E[::m], MAX_DEGREE)
            assert len(actual) == MAX_DEGREE + 1
            if m % 2:
                product_coefficients = convolution(kernel(p, q, m), phi, m * MAX_DEGREE)
                assert actual == [product_coefficients[m * k] for k in range(MAX_DEGREE + 1)]
                odd_identities += 1
            positive = all(min(x) >= 0 for x in actual)
            if premise or m % 2 == 0:
                assert positive, (name, m, actual)
            outcomes.append({'m': m, 'effective_through_degree_8': positive,
                             'phi_coefficients': actual})
        reading = {'name': name, 'dimension': n, 'vertices': vertices,
                   'A': A, 'b': b, 'p': p, 'q': q,
                   'base_effective_through_degree_56': premise,
                   'counts_0_to_8': counts[:9], 'dilates': outcomes}
        all_readings.append(reading)
        print('POLYTOPE', json.dumps(reading, separators=(',', ':')))
    print('SUMMARY', json.dumps({'polytopes': len(all_readings),
          'checked_effective_bases_through_degree_56': effective_premises,
          'ineffective_controls': controls, 'dilation_cases': len(all_readings) * MAX_M,
          'odd_veronese_identities': odd_identities, 'independent_brute_crosschecks': brute_crosschecks}))
    print('BOUNDARY: finite coefficient windows only; no general proof or Lean settlement claimed')
    print('check.py sha256', hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
