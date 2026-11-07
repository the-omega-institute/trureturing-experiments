#!/usr/bin/env python3
"""Exact checks for the Lidbetter--Lin single-look regret conjecture.

No floating point, numerical LP solver, or search-seat encoding is used.
The payoff is reconstructed from Sections 1 and 4 of arXiv:1710.05332v1:
stop upon finding the last object; regret charges only opened empty boxes.
"""
from fractions import Fraction as F
from itertools import combinations, permutations
from math import prod
import sympy as sp


def payoff(c, order, hidden):
    """Return (total search cost, additive regret), by direct search simulation."""
    remaining = set(hidden)
    total = F(0)
    regret = F(0)
    for box in order:
        if not remaining:
            break
        total += c[box]
        if box in remaining:
            remaining.remove(box)
        else:
            regret += c[box]
    assert not remaining
    assert total - sum((c[i] for i in hidden), F(0)) == regret
    return total, regret


def via_positions(c, order, hidden):
    """Independent formulation: empty boxes before the last occupied position."""
    position = {box: j for j, box in enumerate(order)}
    last = max(position[i] for i in hidden)
    return sum((c[i] for i in range(len(c))
                if i not in hidden and position[i] < last), F(0))


def game(c, k):
    placements = tuple(combinations(range(len(c)), k))
    orders = tuple(permutations(range(len(c))))
    matrix = []
    total_matrix = []
    for order in orders:
        entries = [payoff(c, order, h) for h in placements]
        row = [r for t, r in entries]
        assert row == [via_positions(c, order, h) for h in placements]
        matrix.append(row)
        total_matrix.append([t for t, r in entries])
    return placements, orders, matrix, total_matrix


def dot(a, b):
    return sum((x * y for x, y in zip(a, b)), F(0))


def label(items):
    return ''.join(str(i + 1) for i in items)


def rational(x):
    return F(int(x.p), int(x.q))


def main():
    print('MODEL: n! fixed permutations; payoff = additive regret')
    print('ARITHMETIC: Python integers/Fraction and SymPy exact matrices')
    print('SYMPY_VERSION:', sp.__version__)

    # Fidelity: the source's lemma with TeX label le:Wnk, W = k S_(k+1)/S_k.
    for c, k in [((2, 3, 5), 1), ((2, 3, 5), 2), ((19, 13, 13, 4, 2), 2)]:
        h, orders, a, _ = game(c, k)
        weights = [prod(c[i] for i in x) for x in h]
        sk = sum(weights)
        sk1 = sum(prod(c[i] for i in x) for x in combinations(range(len(c)), k + 1))
        expected = F(k * sk1, sk)
        assert all(dot(row, weights) / sk == expected for row in a)
        print('FIDELITY W: costs=%s k=%s permutations=%s S_k=%s W=%s'
              % (c, k, len(orders), sk, expected))

    c = (19, 13, 13, 4, 2)
    h, orders, a, total_a = game(c, 2)
    assert len(h) == 10 and len(orders) == 120
    p = tuple(F(x, 289) for x in (0, 0, 20, 99, 0, 20, 65, 20, 65, 0))
    v = F(4975, 289)
    assert all(x >= 0 for x in p) and sum(p) == 1
    print('PLACEMENT_ORDER:', ' '.join(label(x) for x in h))
    print('HIDER:', ' '.join(map(str, p)))
    print('VALUE:', v)
    for order, row in zip(orders, a):
        lower = dot(row, p)
        assert lower >= v
        print('LOWER', label(order), lower, 'slack', lower - v)
    assert min(dot(row, p) for row in a) == v

    # Independently obtained rational dual, using the eight distinct active rows.
    q = {
        (1, 5, 2, 4, 3): F(128, 11271),
        (1, 5, 3, 4, 2): F(128, 11271),
        (1, 5, 4, 2, 3): F(18076, 146523),
        (1, 5, 4, 3, 2): F(27361, 97682),
        (2, 5, 3, 4, 1): F(14828, 634933),
        (2, 5, 4, 3, 1): F(93133, 293046),
        (3, 5, 2, 4, 1): F(136346, 1904799),
        (3, 5, 4, 2, 1): F(7867, 48841),
    }
    assert len(q) == 8 and sum(q.values()) == 1 and min(q.values()) > 0
    q_rows = [orders.index(tuple(i - 1 for i in order)) for order in q]
    q_weights = list(q.values())
    for order, weight in q.items():
        print('SEARCHER', ''.join(map(str, order)), weight)
    upper = [dot([a[i][j] for i in q_rows], q_weights) for j in range(10)]
    assert all(x <= v for x in upper)
    support = [j for j, x in enumerate(p) if x > 0]
    zero = [j for j, x in enumerate(p) if x == 0]
    assert [j for j, x in enumerate(upper) if x < v] == zero
    for j, x in enumerate(upper):
        print('UPPER', label(h[j]), x, 'slack', v - x)

    # Any optimum is zero on all strictly slack dual columns, and saturates
    # each positive dual row. Normalization plus these equations has rank 6.
    equations = sp.Matrix([[1] * len(support)] +
                          [[a[i][j] for j in support] for i in q_rows])
    rhs = sp.Matrix([1] + [sp.Rational(v.numerator, v.denominator)] * 8)
    assert equations.rank() == len(support) == 6
    solution = equations.gauss_jordan_solve(rhs)[0]
    assert [rational(x) for x in solution] == [p[j] for j in support]
    independent = list(equations.T.rref()[1])
    square = equations[independent, :]
    assert square.det() != 0
    print('UNIQUENESS: zero columns', ' '.join(label(h[j]) for j in zero))
    print('UNIQUENESS: independent equation indices', independent,
          '(0 = normalization; 1..8 = SEARCHER rows)')
    print('UNIQUENESS: matrix', square.tolist())
    print('UNIQUENESS: determinant', square.det(), 'rank', equations.rank())
    assert p[h.index((0, 3))] == p[h.index((1, 3))] > 0
    assert F(c[0], c[1]) == F(19, 13) != 1
    print('RATIO: p14/p24=1; product weighting requires 19/13')

    # Check EVERY nonempty support for the product-weighted strategy.
    # The support is arbitrary, not restricted to a threshold ordering.
    weights = [prod(c[i] for i in hidden) for hidden in h]
    best = None
    best_mask = None
    for mask in range(1, 1 << len(h)):
        w = [weights[j] if (mask >> j) & 1 else 0 for j in range(len(h))]
        z = sum(w)
        assert z > 0
        guaranteed = min(F(dot(row, w), z) for row in a)
        assert guaranteed < v
        if best is None or guaranteed > best:
            best, best_mask = guaranteed, mask
    assert best == F(843, 49) and best_mask == 840
    delta = v - best
    assert delta == F(148, 14161)
    print('ALL_PRODUCT_SUPPORTS: 1023 checked; best', best,
          'support', ' '.join(label(h[j]) for j in range(10) if (best_mask >> j) & 1))
    print('PRODUCT_GAP:', delta)

    # Uniform perturbation certificate, for every real vector with |c'-c|<=r.
    # If c'_i=c_i+x_i, |x_i|<=r<=1, each product changes by <=39r:
    # c'_i*c'_j-c_i*c_j = c_i*x_j+c_j*x_i+x_i*x_j.
    x, y, ci, cj = sp.symbols('x y ci cj')
    assert sp.expand((ci + x) * (cj + y) - ci * cj - (ci * y + cj * x + x * y)) == 0
    r = F(1, 95726048)
    assert 0 < r <= 1 and min(c) - r > 0
    assert max(c) <= 19 and min(weights) == 8 and max(weights) <= 19 * 19
    assert all(len(set(range(5)) - set(hidden)) == 3 for hidden in h)
    assert max(max(row) for row in a) == 45
    # For a support S, W_S>=8, sum |w'_H-w_H|<=10*39r.
    # Normalization gives ||p'_S-p_S||_1 <= 2*10*39r/8 = 195r/2.
    l1_coefficient = F(2 * 10 * 39, 8)
    assert l1_coefficient == F(195, 2)
    product_payoff_coefficient = 3 + 45 * l1_coefficient
    gap_loss_coefficient = 3 + product_payoff_coefficient
    assert gap_loss_coefficient == F(8787, 2)
    margin = delta - gap_loss_coefficient * r
    assert margin == F(97614109, 9381152704) > 0
    print('FAMILY: radius', r, 'l1_coefficient', l1_coefficient,
          'gap_loss_coefficient', gap_loss_coefficient, 'positive_margin', margin)
    print('FAMILY: all real costs in the closed sup-norm cube; all positive rescalings')

    # Expose the payoff discrepancy in the compact candidate, rather than
    # silently treating a regret certificate as a search-cost certificate.
    total_lower = min(dot(row, p) for row in total_a)
    baseline = dot([sum(c[i] for i in hidden) for hidden in h], p)
    assert baseline == F(5169, 289)
    assert total_lower == v + baseline == F(10144, 289)
    assert total_lower != v
    print('TOTAL_COST_DISCREPANCY: occupied-cost expectation', baseline,
          'minimum total cost against candidate', total_lower)
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
