"""Exact checks for Shi–Yau–Zuo, Conjecture 4.79 (DOI 10.4213/im9488e)."""

from itertools import product


def check_binary_ring():
    """Enumerate the full square-zero algebra and generate its ideals by closure."""
    elements = tuple(product(range(2), repeat=4))
    zero = elements[0]
    maximal_ideal = frozenset(a for a in elements if a[0] == 0)

    def add(a, b):
        return tuple(x ^ y for x, y in zip(a, b))

    def multiply(a, b):
        return (a[0] * b[0],) + tuple(
            (a[0] * y) ^ (b[0] * x) for x, y in zip(a[1:], b[1:])
        )

    def ideal(generators):
        values = {zero}
        for g in generators:
            for a in elements:
                scaled = multiply(a, g)
                values |= {add(x, scaled) for x in values}
        return frozenset(values)

    def components(a):
        c, x, y, _ = a
        return (c, x, 0, 0), (c, 0, y, 0), (c, 0, 0, x ^ y)

    def image(values, subset):
        return ideal(g for a in subset for g in values[a])

    def delta(values, target):
        return frozenset(a for a in elements if values[a] <= target)

    identity_values = {a: ideal((a,)) for a in elements}
    values = {a: ideal(components(a)) for a in elements}
    for a, f in product(elements, repeat=2):
        assert identity_values[multiply(a, f)] <= identity_values[f]
        assert values[multiply(a, f)] <= values[f]
    print("Fidelity Example 4.5(1): 256 multiplication-containment tests")

    for images in (identity_values, values):
        targets = set(images.values())
        for target in targets:
            assert image(images, delta(images, target)) == target
        print(f"Fidelity Proposition 4.67: {len(targets)} distinct T-principal ideals are T-full")

    assert len(elements) == 16
    assert len(maximal_ideal) == 8
    assert delta(values, maximal_ideal) == maximal_ideal
    assert image(values, maximal_ideal) == maximal_ideal
    observed = []
    for x, y in product(range(2), repeat=2):
        sizes = {len(values[0, x, y, z]) for z in range(2)}
        assert len(sizes) == 1
        observed.append(sizes.pop())
    assert observed == [1, 4, 4, 4]
    assert all(values[f] != maximal_ideal for f in elements)
    print(f"Binary model: |A|=16, |m|=8, nonunit image sizes={observed}")


def check_prime_field_family(p):
    """Compute in F_p[U,V] using exponent pairs and exact modular coefficients."""

    def multiply(left, right):
        coefficients = {}
        for (i, j), a in left.items():
            for (k, ell), b in right.items():
                degree = i + k, j + ell
                coefficients[degree] = (coefficients.get(degree, 0) + a * b) % p
        return {degree: a for degree, a in coefficients.items() if a}

    def polynomial_product(factors):
        result = {(0, 0): 1}
        for factor in factors:
            result = multiply(result, factor)
        return result

    def evaluate(poly, u, v):
        return sum(a * pow(u, i, p) * pow(v, j, p) for (i, j), a in poly.items()) % p

    finite_factors = [{(1, 0): 1}] + [
        {(1, 0): 1, (0, 1): p - a} for a in range(1, p)
    ]
    factors = [{(0, 1): 1}, *finite_factors]
    witness = polynomial_product(factors)
    assert witness == {(p, 1): 1, (1, p): p - 1}
    for i, factor in enumerate(factors):
        quotient = polynomial_product(factors[:i] + factors[i + 1 :])
        assert multiply(factor, quotient) == witness
    for u, v in product(range(p), repeat=2):
        zero_slots = [i for i, f in enumerate(factors) if evaluate(f, u, v) == 0]
        assert len(zero_slots) == (p + 1 if (u, v) == (0, 0) else 1)
        assert evaluate(witness, u, v) == 0
    assert all(any(evaluate(f, u, v) for u, v in [(1, 0), (0, 1)]) for f in factors)
    print(
        f"F_{p}: {p + 1} coordinate divisibilities, "
        f"{p ** 2} rational pairs checked, nonzero formal polynomial U^{p}V-UV^{p}"
    )


def main():
    check_binary_ring()
    for p in (2, 3, 5, 7):
        check_prime_field_family(p)
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
