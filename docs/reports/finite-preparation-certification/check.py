#!/usr/bin/env python3
"""Exact finite checks for arXiv:2504.02944v3, Conjecture 1.

Reproduce Appendix E, including its warning about theory-dependent witnesses.
Then test finite certificates and response compression on rational qubit data.
These finite tests do not prove the universally quantified conjecture.
Only Python's standard library and rational arithmetic are used.
"""

from fractions import Fraction as Q
from itertools import product


def matrix(rows):
    return tuple(tuple(Q(x) for x in row) for row in rows)


I = matrix([[1, 0], [0, 1]])
X = matrix([[0, 1], [1, 0]])
Z = matrix([[1, 0], [0, -1]])
ZERO = matrix([[0, 0], [0, 0]])


def add(a, b):
    return tuple(tuple(x + y for x, y in zip(r, s)) for r, s in zip(a, b))


def scale(c, a):
    return tuple(tuple(c * x for x in r) for r in a)


def born(rho, effect):
    return sum(rho[i][j] * effect[j][i] for i in range(2) for j in range(2))


def psd(a):
    return a[0][1] == a[1][0] and a[0][0] >= 0 and a[1][1] >= 0 and (
        a[0][0] * a[1][1] - a[0][1] * a[1][0] >= 0
    )


def trace(a):
    return a[0][0] + a[1][1]


def rref(rows):
    a = [list(map(Q, row)) for row in rows]
    pivots = []
    row = 0
    for col in range(len(a[0])):
        pivot = next((j for j in range(row, len(a)) if a[j][col]), None)
        if pivot is None:
            continue
        a[row], a[pivot] = a[pivot], a[row]
        divisor = a[row][col]
        a[row] = [x / divisor for x in a[row]]
        for j in range(len(a)):
            if j != row:
                multiplier = a[j][col]
                a[j] = [x - multiplier * y for x, y in zip(a[j], a[row])]
        pivots.append(col)
        row += 1
        if row == len(a):
            break
    return a, pivots


def identities(operators):
    rows = [[op[i][j] for op in operators] for i in range(2) for j in range(2)]
    a, pivots = rref(rows)
    result = []
    for free in set(range(len(operators))) - set(pivots):
        v = [Q(0)] * len(operators)
        v[free] = Q(1)
        for row, pivot in enumerate(pivots):
            v[pivot] = -a[row][free]
        result.append(v)
    return result


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def check_model(rhos, effects, mu, responses):
    assert all(psd(r) and trace(r) == 1 for r in rhos)
    assert all(psd(e) for e in effects)
    total = ZERO
    for e in effects:
        total = add(total, e)
    assert total == I
    assert all(all(x >= 0 for x in row) and sum(row) == 1 for row in mu)
    assert all(all(x >= 0 for x in r) and sum(r) == 1 for r in responses)
    for rel in identities(rhos):
        assert all(dot(rel, [row[l] for row in mu]) == 0 for l in range(len(responses)))
    for rel in identities(effects):
        assert all(dot(rel, r) == 0 for r in responses)
    for a, rho in enumerate(rhos):
        for b, effect in enumerate(effects):
            assert born(rho, effect) == sum(mu[a][l] * r[b] for l, r in enumerate(responses))


def appendix_e():
    # Source order b=1,2,3,4: Z+, X-, Z-, X+.
    axes = [Z, scale(-1, X), scale(-1, Z), X]
    effects = [scale(Q(1, 4), add(I, axis)) for axis in axes]
    rhos = [scale(Q(1, 2), add(I, scale(-1, axis))) for axis in axes]
    mu = [list(map(Q, row)) for row in [
        [0, Q(1, 2), Q(1, 2), 0], [Q(1, 2), 0, Q(1, 2), 0],
        [Q(1, 2), 0, 0, Q(1, 2)], [0, Q(1, 2), 0, Q(1, 2)],
    ]]
    responses = [list(map(Q, row)) for row in [
        [Q(1, 2), 0, 0, Q(1, 2)], [0, Q(1, 2), Q(1, 2), 0],
        [0, 0, Q(1, 2), Q(1, 2)], [Q(1, 2), Q(1, 2), 0, 0],
    ]]
    check_model(rhos, effects, mu, responses)
    observed = [[born(r, e) for e in effects] for r in rhos]
    expected = [[0, Q(1, 4), Q(1, 2), Q(1, 4)],
                [Q(1, 4), 0, Q(1, 4), Q(1, 2)],
                [Q(1, 2), Q(1, 4), 0, Q(1, 4)],
                [Q(1, 4), Q(1, 2), Q(1, 4), 0]]
    assert observed == expected
    assert sum(observed[b][b] for b in range(4)) == 0
    # Exact sign check: sqrt(2) < 2 because both are positive and 2 < 2^2.
    assert 0 < 2 < 2**2
    print('FIDELITY: Appendix E probability table, ontological model, all identity bases: OK')
    print('FIDELITY: diagonal witness value = 0 < 2-sqrt(2), despite a finite model: OK')


SIGNS = list(product([-1, 1], repeat=2))
RESPONSES = [[Q(1 + s, 4), Q(1 - s, 4), Q(1 + t, 4), Q(1 - t, 4)]
             for s, t in SIGNS]
WEIGHTS = [[u, -u, v, -v] for u, v in SIGNS]


def score(table):
    return sum(dot(w, row) for w, row in zip(WEIGHTS, table))


def certificate():
    # Extreme rays of z>=0, z_--+z_++=z_-++z_+-, normalized by sum z=1.
    # Every ray has one nonzero entry on each side, with equal values.
    rays = []
    for left, right in product([0, 3], [1, 2]):
        z = [Q(0)] * 4
        z[left] = z[right] = Q(1, 2)
        rays.append(z)
    values = [score([[z[a] * r[b] for b in range(4)] for a in range(4)])
              for z, r in product(rays, RESPONSES)]
    assert len(values) == 16 and min(values) == -Q(1, 2) and max(values) == Q(1, 2)
    # Summing ontic columns gives W <= (1/2) sum_a 1 = 2.
    print('CERTIFICATE: all 4 preparation-cone rays x 4 response vertices give W <= 1/2 per unit mass')


def rectangle_checks():
    grid = sorted({Q(n, d) for d in range(1, 7) for n in range(d + 1)})
    counts = {'model': 0, 'violation': 0, 'boundary': 0, 'zero_noise': 0}
    for a, b in product(grid, repeat=2):
        if a * a + b * b > 1:
            continue
        rhos = [scale(Q(1, 2), add(add(I, scale(u * a, X)), scale(v * b, Z)))
                for u, v in SIGNS]
        for eta in grid:
            effects = [scale(Q(1, 4), add(I, scale(eta, axis)))
                       for axis in [X, scale(-1, X), Z, scale(-1, Z)]]
            assert all(psd(r) and trace(r) == 1 for r in rhos)
            assert all(psd(e) for e in effects)
            table = [[born(r, e) for e in effects] for r in rhos]
            assert score(table) == 2 * eta * (a + b)
            if eta == 0:
                check_model(rhos, effects, [[Q(1)] for _ in SIGNS], [[Q(1, 4)] * 4])
                counts['zero_noise'] += 1
            elif eta * (a + b) <= 1:
                mu = [[(1 + eta * (u * a * s + v * b * t)) / 4 for s, t in SIGNS]
                      for u, v in SIGNS]
                check_model(rhos, effects, mu, RESPONSES)
                counts['model'] += 1
                if eta * (a + b) == 1:
                    counts['boundary'] += 1
            else:
                assert score(table) > 2
                counts['violation'] += 1
    assert counts['model'] and counts['violation'] and counts['boundary'] and counts['zero_noise']
    print('RECTANGLES: all a,b,eta in reduced fractions [0,1] with denominator <= 6 and a^2+b^2 <= 1')
    print('RECTANGLES:', counts, 'total =', sum(counts[k] for k in ['model', 'violation', 'zero_noise']))


def finite_detection():
    # Every sampled nonclassical noisy BB84 measurement has a four-pure-state witness.
    pairs = sorted({(Q(q * q - p * p, q * q + p * p), Q(2 * p * q, q * q + p * p))
                    for q in range(1, 41) for p in range(q + 1)})
    etas = sorted({Q(n, d) for d in range(1, 33) for n in range(d + 1)
                   if 2 * n * n > d * d})
    for eta in etas:
        a, b = next((a, b) for a, b in pairs if eta * (a + b) > 1)
        assert a * a + b * b == 1
        assert 2 * eta * (a + b) > 2
    print('FINITE DETECTION:', len(etas), 'distinct rational eta, denominator <= 32, 2*eta^2>1: four pure preparations each')


def compression():
    grid = [Q(-1), -Q(1, 2), Q(0), Q(1, 2), Q(1)]
    interior = list(product(grid, repeat=2))
    responses = [[(1 + x) / 4, (1 - x) / 4, (1 + z) / 4, (1 - z) / 4]
                 for x, z in interior]
    bary = [[(1 + s * x) * (1 + t * z) / 4 for s, t in SIGNS] for x, z in interior]
    assert all(sum(w) == 1 and min(w) >= 0 for w in bary)
    for r, w in zip(responses, bary):
        assert all(r[b] == sum(w[i] * v[b] for i, v in enumerate(RESPONSES)) for b in range(4))
    mu = [[(1 + Q(3, 20) * u * x + Q(1, 5) * v * z) / 25 for x, z in interior]
          for u, v in SIGNS]
    compressed = [[sum(row[l] * bary[l][i] for l in range(25)) for i in range(4)] for row in mu]
    for old, new in zip(mu, compressed):
        assert sum(old) == sum(new) == 1 and min(new) >= 0
        for b in range(4):
            assert sum(old[l] * responses[l][b] for l in range(25)) == sum(new[i] * RESPONSES[i][b] for i in range(4))
    assert all(compressed[0][i] + compressed[3][i] == compressed[1][i] + compressed[2][i] for i in range(4))
    print('COMPRESSION: 25 rational responses -> 4 fixed vertices; statistics and preparation relation preserved: OK')


if __name__ == '__main__':
    appendix_e()
    certificate()
    rectangle_checks()
    finite_detection()
    compression()
    print('ALL_CHECKS_OK')
