#!/usr/bin/env python3
"""Exact certificates for the source coefficient and a finite-pulse counterexample.

Only integers and fractions are used. The finite-pulse certificate checks every
arithmetic hypothesis of the explicitly stated Picard/Gronwall argument in README.
It does not claim to be a Lean proof or to integrate an ODE symbolically.
"""
from fractions import Fraction as Q


def add(p, q):
    return [sum((v[i] if i < len(v) else Q(0) for v in (p, q)), Q(0))
            for i in range(max(len(p), len(q)))]


def mul(p, q):
    out = [Q(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] += a * b
    return out


def scale(c, p):
    return [c * a for a in p]


def evaluate(p, x):
    out = Q(0)
    for c in reversed(p):
        out = out * x + c
    return out


def integral(p, left, right):
    return sum((c * (right ** (i + 1) - left ** (i + 1)) / (i + 1)
                for i, c in enumerate(p)), Q(0))


def source_fidelity():
    first = add([Q(1, 80)], scale(Q(6, 25), mul([0, 1], [Q(1, 4), -1])))
    third = add([Q(1, 80)], scale(Q(15, 4),
                add([1], scale(-1, mul([Q(-25, 3), Q(40, 3)],
                                       [Q(-25, 3), Q(40, 3)])))))
    mean = (integral(first, Q(0), Q(1, 4))
            + Q(1, 80) * (Q(11, 20) - Q(1, 4))
            + integral(third, Q(11, 20), Q(7, 10))
            - 375 * (1 - Q(7, 10)))
    assert mean == Q(-35877, 320)
    assert evaluate(first, Q(1, 8)) == Q(13, 800)
    assert evaluate(third, Q(5, 8)) == Q(301, 80)
    # The arcs are nonnegative additions to 1/80 on their respective intervals:
    assert 40 * Q(11, 20) - 25 == -3
    assert 40 * Q(7, 10) - 25 == 3
    assert evaluate(first, Q(0)) == evaluate(first, Q(1, 4)) == Q(1, 80)
    assert evaluate(third, Q(11, 20)) == evaluate(third, Q(7, 10)) == Q(1, 80)
    print(f'SOURCE integral={mean}; maxima=13/800,301/80; sign switch=7/10')


def kick(state, mass):
    u, v = state
    return u, v - mass * u * u


def drift(state, length):
    u, v = state
    return u + length * v, v


def finite_pulse_certificate():
    b, theta, eta = Q(10000), Q(7, 8), Q(7, 16)
    r, C = 1 / b, b / theta ** 2
    z0, zm, zs, zp = Q(7, 1600), Q(1, 8), Q(1, 2), Q(7, 8)
    assert z0 == eta / 100 and Q(1, 100) ** 2 == r
    points = (z0, zm, zs, zp)
    assert all(0 < x < y < 1 for x, y in zip(points, points[1:]))
    print(f'FAMILY gamma=2; b={b}; r={r}; theta={theta}; C={C}')
    residuals = []
    for z in points:
        s = 1 - z
        D = z - (1 - z) * (z * z + r)
        state = (s, Q(0))
        states = [state]
        for f, x in ((kick, 1), (drift, 1), (kick, b), (drift, r), (kick, -C)):
            state = f(state, x)
            states.append(state)
        assert state[0] == s * D
        F = s ** 2 * (C * D ** 2 - (1 + b * z ** 2))
        assert state[1] == F
        assert all(u > 0 for u, _ in states)
        residuals.append(F)
        print(f'IDEAL z={z}; D={D}; final_v={F}; sign={"+" if F > 0 else "-"}')
    assert [v > 0 for v in residuals] == [False, True, False, True]

    # Polynomial identity proves the uniform lower bound, not merely a grid check:
    # D(z) - (3*z/4-r) = z*(z-1/2)^2+r*z >= 0 for z >= 0.
    Dpoly = [-r, 1 + r, -1, 1]
    lower = [-r, Q(3, 4)]
    square = mul([Q(-1, 2), 1], [Q(-1, 2), 1])
    rhs = add(mul([0, 1], square), [0, r])
    assert add(Dpoly, scale(-1, lower)) == rhs
    delta_D = Q(3, 4) * z0 - r
    delta = (1 - zp) * delta_D
    assert delta_D > 0 and delta > 0
    print(f'UNIFORM ideal_u>={delta}; D>={delta_D}')

    # All five stages are parametrized by tau in [0,1]. On
    # delta/2 <= u <= 2, |v| <= V, their actual vector fields have
    # sup-norm Lipschitz constant <= L, and differ from the ideal
    # fields by at most V*epsilon. See README for the exact fields.
    V, L, exponent = 30000, 60000, 650000
    eps = Q(1, 1 << exponent)
    assert 0 < eps < r < 1
    assert C > 0 and 4 * (C + 1) < L
    assert 1 + b + C < V - 1
    assert eps < 1 and 4 * eps ** 2 <= V * eps
    assert 4 * eps <= V * eps
    B = 1 << (2 * L)  # e^L < 4^L, using the elementary bound e < 4.
    error = Q(0)
    for stage in range(5):
        error = (error + V * eps) * B
        assert error < delta / 2 and error < 1
    assert all(error < abs(f) / 2 for f in residuals)
    assert error <= Q(5 * V, 1 << 50000)
    lengths = (eps, 1 - eps, eps, r - eps, eps)
    coeffs = (1 / eps + eps, eps, b / eps + eps, -eps, -C / eps - eps)
    assert all(t > 0 for t in lengths)
    assert all(c > 0 for c in coeffs[:3])
    assert all(c < 0 for c in coeffs[3:])
    T = sum(lengths, Q(0))
    assert T == 1 + r + eps
    mean = sum((t * c for t, c in zip(lengths, coeffs)), Q(0))
    assert mean == 1 + b - C + eps * (1 + eps - r)
    assert mean < 0
    print(f'ACTUAL epsilon=2^(-{exponent}); all five lengths positive; sign pattern=+++--')
    print(f'GRONWALL L={L}; mismatch<={V}*epsilon; five-stage error<=150000*2^(-50000)')
    print('CERTIFIED margins: positivity, bounded continuation, four strict signs, negative integral')
    print('CONSEQUENCE via local existence, Gronwall, continuous dependence and IVT: >=3 positive Neumann solutions')
    print('SCOPE: exact rational analytic certificates; ODE/IVT bridges are mathematical, not Lean-verified')


if __name__ == '__main__':
    source_fidelity()
    finite_pulse_certificate()
    print('ALL_CHECKS_OK')
