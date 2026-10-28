"""Exact finite checks for PRS arXiv:2610.10745v1, Conjecture 2.

All coefficients are integers or sympy Rational values.  The script checks the
ambient traceful intertwining and the projected traceless identity on several
small (m,n,N) cases.  No floating point arithmetic is used.
"""
from itertools import permutations, product
from fractions import Fraction
import sympy as sp


def compose(p, q):
    """Source convention: (p q)(i) = q(p(i))."""
    return tuple(q[p[i]] for i in range(len(p)))


def inverse(p):
    r = [0] * len(p)
    for i, j in enumerate(p):
        r[j] = i
    return tuple(r)


def transposition(length, i, j):
    p = list(range(length))
    p[i], p[j] = p[j], p[i]
    return tuple(p)


def conjugate(h, g):
    return compose(compose(h, g), inverse(h))


def subgroup(m, n):
    return [tuple(x) + tuple(y[i] + m for i in range(n))
            for x in permutations(range(m)) for y in permutations(range(n))]


def orbit_data(m, n):
    length = m + n
    group = subgroup(m, n)
    remaining = set(permutations(range(length)))
    result = []
    while remaining:
        representative = min(remaining)
        orbit = frozenset(conjugate(h, representative) for h in group)
        result.append((representative, orbit))
        remaining -= orbit
    return result


def variables(size):
    z = sp.Matrix(size, size, lambda i, j: sp.Symbol(f"z{i}{j}"))
    w = sp.Matrix(size, size, lambda i, j: sp.Symbol(f"w{i}{j}"))
    return z, w, list(z) + list(w)


def contraction(g, m, n, size, z, w, traceless=False):
    if traceless:
        z_trace = sum(z[i, i] for i in range(size))
        w_trace = sum(w[i, i] for i in range(size))
        z = sp.Matrix(size, size, lambda i, j:
                      z[i, j] - (z_trace / size if i == j else 0))
        w = sp.Matrix(size, size, lambda i, j:
                      w[i, j] - (w_trace / size if i == j else 0))
    value = 0
    for indices in product(range(size), repeat=m + n):
        term = 1
        for p in range(m):
            term *= z[indices[p], indices[g[p]]]
        for p in range(m, m + n):
            term *= w[indices[p], indices[g[p]]]
        value += term
    return sp.expand(value)


def polynomial_vector(expression, variables_):
    return dict(sp.Poly(expression, *variables_, domain=sp.QQ).terms())


def matrix_of(expressions, variables_, keys=None):
    vectors = [polynomial_vector(e, variables_) for e in expressions]
    if keys is None:
        keys = sorted(set().union(*(set(v) for v in vectors)))
    return sp.Matrix([[v.get(k, 0) for v in vectors] for k in keys]), keys


def h_bar(g, length, size):
    """Definition 2.24 at ξ = 1/N, with source multiplication convention."""
    terms = []
    for i in range(length):
        if i != g[i]:
            t = transposition(length, i, g[i])
            terms += [(Fraction(1), g),
                      (Fraction(-1), compose(compose(t, g), t))]
    for i in range(length):
        for j in range(length):
            if g[j] != i and j != i:
                t = transposition(length, i, g[j])
                u = transposition(length, i, j)
                terms += [(Fraction(1, size), compose(g, t)),
                          (Fraction(-1, size),
                           compose(compose(compose(u, g), t), u))]
    return terms


def permutation_operator(m, n, size, orbits):
    """H on distinct orbit sums, obtained by summing the covariant action."""
    orbit_index = {p: k for k, (_, orbit) in enumerate(orbits) for p in orbit}
    operator = sp.zeros(len(orbits), len(orbits))
    for column, (_, source_orbit) in enumerate(orbits):
        aggregate = {}
        for source in source_orbit:
            for coefficient, target in h_bar(source, m + n, size):
                aggregate[target] = aggregate.get(target, sp.Rational(0)) + sp.Rational(
                    coefficient.numerator, coefficient.denominator)
        for target, coefficient in aggregate.items():
            operator[orbit_index[target], column] = coefficient
    return operator


def delta_prime(expression, size, z, w):
    z_trace = sum(z[i, i] for i in range(size))
    w_trace = sum(w[i, i] for i in range(size))
    z_prime = sp.Matrix(size, size, lambda i, j:
                        z[i, j] - (z_trace / size if i == j else 0))
    w_prime = sp.Matrix(size, size, lambda i, j:
                        w[i, j] - (w_trace / size if i == j else 0))

    def dz(i, j, e):
        return sp.diff(e, z[i, j]) - (sum(sp.diff(e, z[k, k]) for k in range(size)) / size
                                      if i == j else 0)

    def dw(i, j, e):
        return sp.diff(e, w[i, j]) - (sum(sp.diff(e, w[k, k]) for k in range(size)) / size
                                      if i == j else 0)

    commutator = w_prime * z_prime - z_prime * w_prime
    result = 0
    for i in range(size):
        for j in range(size):
            derivative_commutator = 0
            for k in range(size):
                derivative_commutator += dz(i, k, dw(k, j, expression))
                derivative_commutator -= dw(i, k, dz(k, j, expression))
            result += commutator[i, j] * derivative_commutator
    return sp.expand(-result / size)


def check_case(m, n, size):
    orbits = orbit_data(m, n)
    z, w, variables_ = variables(size)
    orbit_contractions = [len(orbit) * contraction(rep, m, n, size, z, w)
                          for rep, orbit in orbits]
    traceless_contractions = [len(orbit) * contraction(rep, m, n, size, z, w, True)
                              for rep, orbit in orbits]
    q_matrix, _ = matrix_of(orbit_contractions, variables_)
    metric = sp.diag(*[len(orbit) for _, orbit in orbits])
    kernel = (sp.Matrix.hstack(*q_matrix.nullspace())
              if q_matrix.nullspace() else sp.zeros(len(orbits), 0))

    derangement_columns = [k for k, (rep, _) in enumerate(orbits)
                           if all(rep[i] != i for i in range(m + n))]
    constraints = kernel.T * metric[:, derangement_columns]
    intersection_coordinates = (sp.Matrix.hstack(*constraints.nullspace())
                                 if constraints.nullspace()
                                 else sp.zeros(len(derangement_columns), 0))
    intersection = sp.zeros(len(orbits), intersection_coordinates.cols)
    for row, column in enumerate(derangement_columns):
        intersection[column, :] = intersection_coordinates[row, :]
    projection = (intersection * (intersection.T * metric * intersection).inv()
                  * intersection.T * metric
                  if intersection.cols else sp.zeros(len(orbits), len(orbits)))
    operator = permutation_operator(m, n, size, orbits)

    # Include every monomial from both sides so no output term can be dropped.
    traceless_deltas = [delta_prime(e, size, z, w) for e in traceless_contractions]
    keys = sorted(set().union(*(set(polynomial_vector(e, variables_))
                                for e in traceless_contractions + traceless_deltas)))
    f_matrix, _ = matrix_of(traceless_contractions, variables_, keys)
    delta_matrix, _ = matrix_of(traceless_deltas, variables_, keys)

    # Fidelity: q H = Δ q on the complete PCA.  Matrix derivatives in
    # the paper are transposed: (∂Z)_{ji} = ∂/∂Z_{ij}.  Transposing the
    # derivative commutator reverses its order; the following coordinate
    # formula is exactly -tr([Z,W][∂Z,∂W])/N, with normal ordering.
    def traceful_delta(e):
        commutator = w * z - z * w
        result = 0
        for i in range(size):
            for j in range(size):
                derivative_commutator = 0
                for k in range(size):
                    derivative_commutator += sp.diff(sp.diff(e, w[k, j]), z[i, k])
                    derivative_commutator -= sp.diff(sp.diff(e, z[k, j]), w[i, k])
                result += commutator[i, j] * derivative_commutator
        return sp.expand(-result / size)
    traceful_deltas = [traceful_delta(e) for e in orbit_contractions]
    q_keys = sorted(set().union(*(set(polynomial_vector(e, variables_))
                                 for e in orbit_contractions + traceful_deltas)))
    q_mat, _ = matrix_of(orbit_contractions, variables_, q_keys)
    traceful_delta_matrix, _ = matrix_of(traceful_deltas, variables_, q_keys)
    fidelity_ok = q_mat * operator == traceful_delta_matrix
    assert f_matrix * projection == f_matrix, "fP = f failed"

    projected_ok = f_matrix * projection * operator * intersection == delta_matrix * intersection
    return fidelity_ok, projected_ok, intersection.cols, len(orbits), q_matrix.rank()


def main():
    fidelity, projected, dimension, orbit_count, q_rank = check_case(2, 2, 2)
    assert fidelity, "traceful source intertwining failed"
    print(f"fidelity (m,n,N)=(2,2,2): ok; orbit_count={orbit_count}; q_rank={q_rank}")
    cases = [(0, 0, 1), (0, 0, 2), (1, 0, 1), (0, 1, 2),
             (1, 1, 1), (1, 1, 2), (2, 0, 2), (0, 2, 2),
             (2, 1, 2), (2, 2, 2), (2, 2, 3), (3, 1, 2)]
    for case in cases:
        result = check_case(*case)
        assert result[0] and result[1], f"identity failed for {case}: {result}"
        print(f"conditional identity {case}: ok; U_dimension={result[2]}")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
