"""Exact arithmetic check for the HKMP positive-weight K4 family."""

from itertools import combinations
from sympy import Matrix, Rational, Symbol, cancel, expand


t = Symbol("t", real=True)
x, y = Symbol("x", real=True), Symbol("y", real=True)
vertices = range(4)
edge_order = list(combinations(vertices, 2))
weights = dict(zip(edge_order, (1, 2, 3, 3, t, 2)))


def weight(i, j):
    return weights[tuple(sorted((i, j)))]


def neumann_pair_from_source_form(i, j):
    """Evaluate the source's Q_N / ||f||^2 on the mean-zero pair (x,-x)."""
    value = weight(i, j) * (2 * x) ** 2
    for z in set(vertices) - {i, j}:
        a, b = weight(i, z), weight(j, z)
        harmonic = (a * x - b * x) / (a + b)
        value += a * (x - harmonic) ** 2 + b * (-x - harmonic) ** 2
    return cancel(value / (2 * x * x))


expected_neumann = {
    (0, 1): Rational(22, 5) + 6 * t / (t + 3),
    (0, 2): Rational(79, 10),
    (0, 3): 8 + 2 * t / (t + 1),
    (1, 2): Rational(22, 3) + 4 * t / (t + 2),
    (1, 3): 2 * t + Rational(39, 10),
    (2, 3): Rational(32, 5) + 6 * t / (t + 3),
}

for pair, expected in expected_neumann.items():
    assert cancel(neumann_pair_from_source_form(*pair) - expected) == 0

# An unordered exhaustive two-partition is represented by the side containing 0.
partitions = []
for mask in range(1, 15):
    left = frozenset(i for i in vertices if mask & (1 << i))
    right = frozenset(vertices) - left
    if 0 in left:
        partitions.append((left, right))
assert len(partitions) == 7
assert sorted(len(left) for left, _ in partitions) == [1, 2, 2, 2, 3, 3, 3]
assert sum(min(len(left), len(right)) == 1 for left, right in partitions) == 4

# The three pairings in the source argument are A={(0,1),(2,3)},
# C={(0,2),(1,3)}, and P={(0,3),(1,2)}.
N_A = expected_neumann[(2, 3)]
N_C = expected_neumann[(1, 3)]
N_P = expected_neumann[(1, 2)]
assert cancel((N_A - N_P) - (16 * (t - 6) ** 2 + 122 * (t - 6) + 72) /
              (15 * (t + 2) * (t + 3))) == 0
assert cancel(N_P - Rational(31, 3) - (t - 6) / (t + 2)) == 0
assert cancel(10 - expected_neumann[(0, 3)] - 2 / (t + 1)) == 0
assert cancel(N_C - Rational(159, 10) - 2 * (t - 6)) == 0
assert cancel(expected_neumann[(2, 3)] - expected_neumann[(0, 1)]) == 2
assert cancel(expected_neumann[(1, 3)] - expected_neumann[(0, 2)] - 2 * (t - 2)) == 0

# Boundaryless energies are 2*w_ij for a two-vertex mean-zero cluster.
assert max(2 * weight(0, 3), 2 * weight(1, 2)) == 6
assert max(2 * weight(0, 1), 2 * weight(2, 3)) == 4

# Dirichlet certificates from the source proof.
dirichlet_pair = Matrix([[t + 4, -3], [-3, 7]])
lhs = (t + 4) * x * x - 6 * x * y + 7 * y * y
rhs = Rational(51, 10) * (x * x + y * y) + (t - 6) * x * x + ((19 * y - 30 * x) ** 2 + 31 * x * x) / 190
assert expand(lhs - rhs) == 0
competitor_02 = Matrix([[6, -2], [-2, 7]])
competitor_13 = Matrix([[t + 4, -t], [-t, t + 5]])
u = Matrix([2, 1])
v = Matrix([1, 1])
assert (u.T * competitor_02 * u)[0] / (u.T * u)[0] == Rational(23, 5)
assert cancel((v.T * competitor_13 * v)[0] / (v.T * v)[0] - Rational(9, 2)) == 0
assert Rational(9, 2) < Rational(23, 5) < Rational(51, 10)

print("Exact source-form Neumann pair energies:")
for pair, value in expected_neumann.items():
    print(f"  {pair}: {value.subs(t, 6) if hasattr(value, 'subs') else value}")
print("At t=6, E_N(A), E_N(C), E_N(P):", [N_A.subs(t, 6), N_C.subs(t, 6), N_P.subs(t, 6)])
print("The seven unordered exhaustive bipartitions and all exact comparison certificates passed.")
print("ALL_CHECKS_OK")
