"""Finite Arb certificates for pure and mixed rate tomography on the charged shell."""
from functools import lru_cache, reduce
from itertools import combinations, combinations_with_replacement
from math import prod
from pathlib import Path
from runpy import run_path
import json

from flint import arb, arb_mat

root = Path(__file__).resolve().parent
model = run_path(str(root / 'ch39-shell-optimal-static-readings.py'))
data = json.loads((root / 'ch39-shell-arc-rate-certificate.json').read_text())
pairs = list(combinations_with_replacement(range(7), 2))
spread = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 18, 22, 27, 34,
          44, 58, 78, 107, 148, 210, 315, 470, 730, 1150, 1900, 3200]


@lru_cache(None)
def response(i):
    return model['row'](model['shift'](i) if i else arb(0), True)


def measurement_matrix(indices):
    return arb_mat([[r[i] * r[j] * (1 if i == j else 2) for i, j in pairs]
                    for r in map(response, indices)])


def block(matrix, rows, columns):
    return arb_mat([[matrix[i, j] for j in columns] for i in rows])


def row_sum_norm(matrix):
    return reduce(arb.max, (sum(abs(x) for x in row) for row in matrix.tolist()))


def mixed_design(indices, lower, upper):
    matrix = measurement_matrix(indices)
    inverse = matrix.inv()
    factor = row_sum_norm(inverse)
    assert arb(lower) < factor < arb(upper)
    print('Mixed design', list(indices), 'sharp entry error', factor.str(12),
          'condition number', (factor * row_sum_norm(matrix)).str(12))


def allowed_vectors(equations, independent_rows, pivots):
    """Parametrize the solutions of the selected homogeneous equations."""
    n = equations.ncols()
    free = [j for j in range(n) if j not in pivots]
    pivot_block = block(equations, independent_rows, pivots)
    if not free:
        assert abs(pivot_block.det()) > 0
        return None
    solved = pivot_block.solve(-block(equations, independent_rows, free))
    basis = arb_mat(n, len(free))
    for k, j in enumerate(free):
        basis[j, k] = 1
        for i, pivot in enumerate(pivots):
            basis[pivot, k] = solved[i, k]
    return basis


def kernel_matrices(indices, order):
    # The witnesses name their free variables in the supplied column order.
    assert sorted(order) == list(range(28))
    matrix = measurement_matrix(indices)
    reordered = block(matrix, range(len(indices)), order)
    basis = allowed_vectors(reordered, range(len(indices)), list(range(len(indices))))
    result = []
    for k in range(basis.ncols()):
        form = arb_mat(7, 7)
        for l, pair_id in enumerate(order):
            i, j = pairs[pair_id]
            form[i, j] = form[j, i] = basis[l, k]
        result.append(form)
    return result


def monomials(n, degree):
    if n == 1:
        yield (degree,)
        return
    for i in range(degree + 1):
        for rest in monomials(n - 1, degree - i):
            yield (i,) + rest


def minor_coefficients(forms, powers):
    """Interpolate homogeneous degree-five coefficients at the integer simplex."""
    minors = list(combinations_with_replacement(list(combinations(range(7), 5)), 2))
    evaluations = []
    for point in powers:
        form = sum((a * h for a, h in zip(point, forms, strict=True)), arb_mat(7, 7))
        evaluations.append([block(form, rows, cols).det() for rows, cols in minors])
    vandermonde = arb_mat([[prod(t ** k for t, k in zip(point, alpha, strict=True))
                            for alpha in powers] for point in powers])
    return vandermonde.solve(arb_mat(evaluations)).transpose()


def extend_degree(powers, basis, stage):
    """At a common zero, t_j times its degree-d evaluation lies in this basis."""
    n = len(powers[0])
    dim = basis.ncols()
    following = list(monomials(n, sum(powers[0]) + 1))
    lookup = {a: i for i, a in enumerate(powers)}
    first = []
    equalities = []
    for alpha in following:
        parents = []
        for j, power in enumerate(alpha):
            if power:
                beta = list(alpha)
                beta[j] -= 1
                parents.append((j, lookup[tuple(beta)]))
        first.append(parents[0])
        equalities.extend((parents[0], other) for other in parents[1:])

    # Only the listed equalities are needed. Propagate a containing space.
    equations = arb_mat(len(stage['rows']), n * dim)
    for row, index in enumerate(stage['rows']):
        (i, u), (j, v) = equalities[index]
        for k in range(dim):
            equations[row, i * dim + k] = basis[u, k]
            equations[row, j * dim + k] = -basis[v, k]
    lifted = allowed_vectors(equations, range(equations.nrows()), stage['columns'])
    if lifted is None:
        return following, None
    result = arb_mat([[sum(basis[u, k] * lifted[i * dim + k, l] for k in range(dim))
                       for l in range(lifted.ncols())] for i, u in first])
    return following, result


def pure_design(indices):
    forms = kernel_matrices(indices, data['columns'])
    assert len(forms) == 6
    powers = list(monomials(len(forms), 5))
    stage, *rest = data['elimination']
    basis = allowed_vectors(minor_coefficients(forms, powers), stage['rows'], stage['columns'])
    assert basis.ncols() == 56
    for step, dim in zip(rest, (21, 6, 1, 0), strict=True):
        powers, basis = extend_degree(powers, basis, step)
        assert (0 if basis is None else basis.ncols()) == dim
    assert sum(powers[0]) == 9
    print('Pure design', indices, 'degree 5-9 evaluation space bounds', [56, 21, 6, 1, 0])


def schur_equations(forms, fixed, point):
    free_forms = [form for j, form in enumerate(forms) if j != fixed]
    form = forms[fixed] + sum((value * h for value, h in zip(point, free_forms, strict=True)),
                              arb_mat(7, 7))
    lead = block(form, range(4), range(4))
    cross = block(form, range(4), range(4, 7))
    lower = block(form, range(4, 7), range(4, 7))
    tangent = lead.solve(-cross)
    schur = lower + cross.transpose() * tangent
    stack = arb_mat(tangent.tolist() + [[int(i == j) for j in range(3)] for i in range(3)])
    derivatives = [stack.transpose() * h * stack for h in free_forms]
    entries = list(combinations_with_replacement(range(3), 2))
    values = arb_mat([[schur[i, j]] for i, j in entries])
    jacobian = arb_mat([[d[i, j] for d in derivatives] for i, j in entries])
    return form, lead, values, jacobian


def deletion(indices, witness):
    omitted = witness['omitted']
    remaining = [i for i in indices if i != omitted]
    permutation = witness['permutation']
    assert sorted(permutation) == list(range(7))
    forms = [block(h, permutation, permutation) for h in kernel_matrices(remaining, witness['columns'])]
    center = [arb(v) for v in witness['center']]
    radius = arb('1e-55')
    box = [v + arb(0, radius) for v in center]
    _, _, values, jacobian = schur_equations(forms, witness['fixed'], center)
    form, lead, _, enclosure = schur_equations(forms, witness['fixed'], box)
    preconditioner = jacobian.mid().inv().mid()
    identity = arb_mat([[int(i == j) for j in range(6)] for i in range(6)])
    error = identity - preconditioner * enclosure
    # The map t -> t - Y F(t) contracts this closed box into itself.
    assert row_sum_norm(preconditioner * values) + radius * row_sum_norm(error) < radius
    minors = [arb(1)] + [block(lead, range(i), range(i)).det() for i in range(1, 5)]
    assert all(abs(v) > 0 for v in minors)
    assert sum(a * b < 0 for a, b in zip(minors, minors[1:])) == 2

    gaps = []
    for index in (0, omitted):
        r = response(index)
        column = arb_mat([[r[j]] for j in permutation])
        gap = (column.transpose() * form * column)[0, 0]
        assert abs(gap) > 0
        gaps.append(abs(gap))
    entry_gap = reduce(arb.max, (abs(v) for v in form.entries()))
    return gaps[0] / gaps[1], entry_gap / gaps[1]


def main():
    mixed_design(range(1, 29), '4.65043e100', '4.65045e100')
    mixed_design(spread, '5.97595e73', '5.97596e73')
    indices = data['indices']
    assert len(set(indices)) == len(indices) == 22
    pure_design(indices)
    assert [w['omitted'] for w in data['deletions']] == indices
    lower_bounds = [deletion(indices, w) for w in data['deletions']]
    best_static, best_entry = (reduce(arb.max, col) for col in zip(*lower_bounds, strict=True))
    assert best_static > arb("66110.96") and best_entry > arb("2.0619e43")
    print('All 22 single deletions admit signature (2,2) collisions with a static rate gap.')
    print('Pure design error factors at least', best_static.str(12), '(static rate),',
          best_entry.str(12), '(matrix entries)')


if __name__ == '__main__':
    main()
