"""Exact checks for arXiv:2610.24431v1, §4.1 and equations (18)–(19)."""

from sympy import Matrix, Rational, diff, diag, expand, factor, ones, sign, symbols


def require(condition, description):
    if not condition:
        raise AssertionError(description)


def same_vector(left, right, description):
    require(left.shape == right.shape, description)
    require(all(factor(x - y) == 0 for x, y in zip(left, right)), description)


def local_cost(t):
    return (t - 1) ** 2 / 2


def objective(x):
    return sum(local_cost(value) for value in x)


def signed_gradient(x):
    return x.applyfunc(lambda value: sign(factor(value - 1)))


def step(x, tau, laplacian):
    return x - tau * laplacian * signed_gradient(x)


def main():
    a, delta, t, tau, z0, z1 = symbols("a delta t tau z0 z1", real=True)
    weights = Matrix([[0, 1], [1, 0]])
    one = ones(2, 1)
    degree = weights * one
    laplacian = diag(*degree) - weights
    require(weights.T * one == degree, "weight balance")
    require(all(weights[i, i] == 0 and weights[i, 1 - i] > 0 for i in range(2)),
            "no loops and both directed edges")
    require(laplacian.eigenvals() == {0: 1, 2: 1}, "Laplacian spectrum")
    require(diff(local_cost(t), t) == t - 1, "local gradient")
    require(diff(local_cost(t), t, 2) == 1, "cost curvature")

    # Equation (19) with an arbitrary input and signal, covering every sign case.
    x = Matrix([t, 4 - t])
    signal = Matrix([z0, z1])
    rhs = x - tau * laplacian * signal
    require(expand(sum(rhs) - 4) == 0, "source anytime feasibility")
    same_vector(rhs, Matrix([t - tau * (z0 - z1), 4 - t - tau * (z1 - z0)]),
                "equations (18) and (19) agree")
    same_vector(x - tau * laplacian * one, x, "Lemma 7 stationary sign consensus")
    print("FIDELITY_OK equations (18)-(19), anytime feasibility, Lemma 7")

    x0 = Matrix([1 - delta, 3 + delta])
    x_star = Matrix([2, 2])
    require(expand(objective(x) - objective(x_star) - (t - 2) ** 2) == 0,
            "unique constrained optimizer")
    # Parameterize the whole open triangle 0 < delta < 2*a < 1.
    p, q, r, s = symbols("p q r s", positive=True)
    family = {a: p / (2 * (p + q)), delta: p * r / ((p + q) * (r + s))}
    small_quantities = [
        delta, 1 - delta, 2 * a - delta, 1 - 2 * a + delta, 1 - 2 * a,
    ]
    require(all(factor(value.subs(family)).is_positive is True
                for value in small_quantities), "open family parameter margins")
    require(signed_gradient(x0.subs(family)) == Matrix([-1, 1]),
            "initial sign disagreement")

    u = 2 * a - delta
    x1 = Matrix([1 + u, 3 - u])
    same_vector(step(x0.subs(family), family[a], laplacian), x1.subs(family),
                "first harmonic update")
    require(signed_gradient(x1.subs(family)) == one, "strictly positive final gradients")
    same_vector(step(x1.subs(family), tau, laplacian), x1.subs(family),
                "stationary for an arbitrary later step")
    require(all(factor(value.subs(family)).is_positive is True
                for state in (x0, x1, x_star) for coord in state
                for value in (coord, 4 - coord)), "strictly inside the boxes")
    require(expand(objective(x1) - objective(x_star) - (1 - u) ** 2) == 0,
            "positive objective gap formula")
    print("FAMILY_OK 0 < delta < 2*a < 1, one-step entry, arbitrary later steps")

    sample = {a: Rational(1, 100), delta: Rational(1, 100)}
    start = x0.subs(sample)
    after = step(start, sample[a], laplacian)
    require(after == Matrix([Rational(101, 100), Rational(299, 100)]), "sample state")
    gap = objective(after) - objective(x_star)
    require(gap == Rational(9801, 10000), "sample objective gap")
    distance_squared = (after - x_star).dot(after - x_star)
    require(distance_squared == 2 * gap, "sample squared Euclidean distance")
    for k in range(1, 33):
        require(step(after, sample[a] / (k + 1), laplacian) == after,
                f"sample harmonic step {k + 1}")
    print(f"EXAMPLE x0={tuple(map(str, start))} x1={tuple(map(str, after))}")
    print(f"EXAMPLE cost={objective(after)} optimum={objective(x_star)} gap={gap}")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
