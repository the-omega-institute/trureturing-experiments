"""Exact CP-conjugate witnesses for a normalized three-flavor trace experiment."""
from itertools import product

from sympy import I, Rational, expand, eye, im, ones


def reflections(basis):
    return [
        (eye(3) - 2 * basis[:, j] * basis[:, j].H).applyfunc(expand)
        for j in range(3)
    ]


def verify(sign):
    v = eye(3) + (sign * I - 1) * ones(3) / 3
    assert (v.H * v).applyfunc(expand) == eye(3)
    up, down = reflections(eye(3)), reflections(v)
    moduli = v.applyfunc(lambda z: expand(z * z.conjugate()))
    for i, j in product(range(3), repeat=2):
        assert expand((up[i] * down[j]).trace()) == 4 * moduli[i, j] - 1

    for a, b, c in product(up + down, repeat=3):
        assert im(expand((a * b * c).trace())) == 0

    for i, j, k, l in product(range(3), repeat=4):
        u = up[i] * down[j] * up[k] * down[l]
        quartet = v[i, j] * v[k, j].conjugate() * v[k, l] * v[i, l].conjugate()
        assert im(expand(u.trace())) == 16 * im(expand(quartet))

    u = (up[0] * down[0] * up[1] * down[1]).applyfunc(expand)
    j_cp = sign * Rational(2, 27)
    # The three sides are |V_0j conjugate(V_1j)|, evaluated as squares.
    a, b, c = (moduli[0, j] * moduli[1, j] for j in range(3))
    assert a * b - ((c - a - b) / 2)**2 == j_cp**2
    assert im(expand(u.trace())) / 3 == 16 * j_cp / 3

    # Project the control onto (|0> + i |1>) / sqrt(2), then sum over flavors.
    amplitude = eye(3) - I * u
    prob_plus = expand((amplitude.H * amplitude).trace()) / 12
    assert prob_plus == (1 + 16 * j_cp / 3) / 2
    return prob_plus


if __name__ == "__main__":
    pair = verify(1), verify(-1)
    assert pair == (Rational(113, 162), Rational(49, 162))
    print("CP-conjugate Y+ probabilities:", *pair)
