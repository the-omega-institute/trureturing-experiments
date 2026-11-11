#!/usr/bin/env python3
"""Exact checks for the positive/alternating Lüroth continuity candidate."""

from fractions import Fraction as Q
from itertools import product


def branch(d, alternating=False):
    assert isinstance(d, int) and d >= 2
    scale = Q(1, d * (d - 1))
    return (Q(1, d - 1), -scale) if alternating else (Q(1, d), scale)


def compose(f, g):
    return (f[0] + f[1] * g[0], f[1] * g[1])


def word(ds, alternating=False):
    f = (Q(0), Q(1))
    for d in ds:
        f = compose(f, branch(d, alternating))
    return f


def apply(f, x):
    return f[0] + f[1] * x


def periodic_value(prefix, period, alternating=False):
    cycle = word(period, alternating)
    tail = cycle[0] / (1 - cycle[1])
    return apply(word(prefix, alternating), tail)


def alt_step(x, d):
    # Strict containment in the open alternating branch interval.
    assert Q(1, d) < x < Q(1, d - 1)
    return (Q(1, d - 1) - x) * d * (d - 1)


def source_term(ds, n, alternating=False):
    denominator = 1
    for d in ds[:n]:
        denominator *= d * (d - 1)
    denominator *= ds[n] - 1 if alternating else ds[n]
    return Q((-1) ** n if alternating else 1, denominator)


def main():
    # Source fidelity: Moroz v2, Proposition 2.5, for phi_n = 1.
    count = 0
    for k in range(1, 5):
        for ds in product(range(2, 9), repeat=k):
            pos, alt = word(ds), word(ds, True)
            diameter = Q(1)
            for d in ds:
                diameter /= d * (d - 1)
            assert pos[1] == abs(alt[1]) == diameter
            assert pos[0] == sum(source_term(ds, n) for n in range(k))
            assert alt[0] == sum(source_term(ds, n, True) for n in range(k))
            count += 1
    print(f"FIDELITY Proposition 2.5: {count} words; equal cylinder diameters")
    print("FIDELITY Section 2: affine branches equal the printed series partial sums")

    assert periodic_value([], [2]) == 1
    assert periodic_value([], [2], True) == Q(2, 3)
    b = periodic_value([2, 8], [2])
    fb = periodic_value([2, 8], [2], True)
    right = apply(word([2, 7], True), Q(0))
    assert b == Q(4, 7)
    assert fb == Q(157, 168)
    assert right == Q(11, 12)
    assert fb - right == Q(1, 56)
    x = b
    orbit = [x]
    for d in [2, 2, 4]:
        x = alt_step(x, d)
        orbit.append(x)
    assert orbit == [Q(4, 7), Q(6, 7), Q(2, 7), Q(4, 7)]
    assert periodic_value([], [2, 2, 4], True) == b
    print(f"BASE b={b} F(b)={fb} right_limit={right} jump={fb-right}")
    print("BASE alternating orbit=" + " -> ".join(map(str, orbit)))

    # Exact sequences and range bounds support the proposed limit argument.
    for n in range(2, 101):
        t = periodic_value([n], [2])
        ft = periodic_value([n], [2], True)
        xn = apply(word([2, 7]), t)
        yn = apply(word([2, 7], True), ft)
        assert t == Q(1, n - 1)
        assert xn == b + Q(1, 84 * (n - 1))
        assert 0 < ft < Q(1, n - 1)
        assert 0 < yn - right < Q(1, 84 * (n - 1))
        assert xn > b
    print("BASE exact approach and output error bounds: N=2..100")

    for j in range(2, 21):
        bj = Q(7 * j - 4, 7 * j * (j - 1))
        assert bj == apply(branch(j, True), b)
        assert bj == apply(word([j, 3]), b)
        assert bj == periodic_value([j, 3, 2, 8], [2])
        assert bj == periodic_value([j], [2, 2, 4], True)
        assert alt_step(bj, j) == b
        fj = periodic_value([j, 3, 2, 8], [2], True)
        rj = apply(word([j, 3], True), right)
        jump = Q(1, 336 * j * (j - 1))
        assert fj == apply(word([j, 3], True), fb)
        assert fj - rj == jump > 0
        assert 0 < bj <= Q(1, j - 1) <= 1
        next_b = Q(7 * (j + 1) - 4, 7 * (j + 1) * j)
        assert bj - next_b == Q(7 * j - 1, 7 * j * (j - 1) * (j + 1)) > 0
        for n in range(2, 41):
            xn = periodic_value([j, 3, 2, 7, n], [2])
            yn = periodic_value([j, 3, 2, 7, n], [2], True)
            bound = Q(1, 504 * j * (j - 1) * (n - 1))
            assert xn == bj + bound
            assert 0 < yn - rj < bound
        print(f"FAMILY j={j} b_j={bj} F(b_j)={fj} right_limit={rj} jump={jump} "
              "alternating_code=(j,(2,2,4)^infinity) interior=True")
    print("No floating point arithmetic; limit bounds are exact rational identities.")
    print("These checks are finite evidence, not a Lean proof of any limit or universal claim.")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
