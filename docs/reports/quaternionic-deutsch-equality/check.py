#!/usr/bin/env python3
"""Exact finite checks for Krishna arXiv:2610.35881v1, Theorem 2.1 and Q2.12.

No floats, numerical logarithms, random sampling, or external dependencies.
These finite checks are evidence about the stated families, not a proof of the
classification over arbitrary continuous frames.
"""
from fractions import Fraction as F
from itertools import product
from math import lcm


def quat(a=0, b=0, c=0, d=0):
    return tuple(map(F, (a, b, c, d)))


ZERO, ONE, I, J, K = quat(), quat(1), quat(0, 1), quat(0, 0, 1), quat(0, 0, 0, 1)


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def scale(r, a):
    return tuple(r * x for x in a)


def conj(a):
    return (a[0], -a[1], -a[2], -a[3])


def mul(a, b):
    w, x, y, z = a
    t, u, v, s = b
    return (w*t-x*u-y*v-z*s, w*u+x*t+y*s-z*v,
            w*v-x*s+y*t+z*u, w*s+x*v-y*u+z*t)


def norm2(a):
    return sum(x*x for x in a)


def inner(a, b):
    r = ZERO
    for x, y in zip(a, b):
        r = add(r, mul(conj(x), y))
    return r


def vnorm2(a):
    return sum(norm2(x) for x in a)


def basis(d, rotated):
    # A real rational rotation with quaternionic phases on both sides.
    rows = [ONE, I, J]
    cols = [J, K, I]
    matrix = [[F(i == j) for j in range(d)] for i in range(d)]
    if rotated and d >= 2:
        matrix[0][0], matrix[0][1] = F(3, 5), F(-4, 5)
        matrix[1][0], matrix[1][1] = F(4, 5), F(3, 5)
    return [tuple(scale(matrix[i][j], mul(rows[i], cols[j])) for i in range(d))
            for j in range(d)]


def frame(d, rotated, r, copies):
    b = basis(d, rotated)
    return [(F(1, copies) / (r*r), tuple(scale(r, x) for x in a))
            for a in b for _ in range(copies)]


def parseval_matrix(fr, d):
    # Sum w*a*a^* = Id, so the check is for every quaternionic vector.
    for i in range(d):
        for j in range(d):
            z = ZERO
            for w, a in fr:
                z = add(z, scale(w, mul(a[i], conj(a[j]))))
            assert z == (ONE if i == j else ZERO), (i, j, z)
    assert all(w >= 0 and vnorm2(a) <= 1 for w, a in fr)


def buzano(a, b, x):
    # 2*sqrt(L) <= X*(sqrt(C)+sqrt(A*B)), decided by rational squaring.
    X, A, B = vnorm2(x), vnorm2(a), vnorm2(b)
    L = norm2(inner(a, x)) * norm2(inner(x, b))
    C = norm2(inner(a, b))
    if not X:
        return L == 0
    t = 4*L/(X*X) - C - A*B
    return t <= 0 or t*t <= 4*C*A*B


def qmul(x, y, D):
    a, b = x
    c, d = y
    return (a*c+b*d*D, a*d+b*c)


def qpow(x, n, D):
    z = (F(1), F(0))
    while n:
        if n & 1:
            z = qmul(z, x, D)
        x = qmul(x, x, D)
        n //= 2
    return z


def compare_quad_rational(x, q, D):
    # Sign of a+b*sqrt(D)-q. Our bases and powers have b,D >= 0.
    a, b = x
    assert b >= 0 and D >= 0
    if b == 0 or D == 0:
        return (a > q) - (a < q)
    if a >= q:
        return 1
    t = b*b*D-(q-a)**2
    return (t > 0) - (t < 0)


def state_check(tau, omega, x, d):
    n = vnorm2(x)
    assert n > 0
    ps = [norm2(inner(a, x))/n for _, a in tau]
    qs = [norm2(inner(a, x))/n for _, a in omega]
    if any(p == 0 for p in ps+qs):
        return None
    assert all(0 < p <= 1 for p in ps+qs)
    assert sum(w*p for (w, _), p in zip(tau, ps)) == 1
    assert sum(w*q for (w, _), q in zip(omega, qs)) == 1
    D = max(norm2(inner(a, b)) for _, a in tau for _, b in omega)
    assert 0 <= D <= 1
    # M=((1+c)/2)^2, c=sqrt(D).
    M = ((1+D)/4, F(1, 2))
    A = max(p for (w, _), p in zip(tau, ps) if w)
    B = max(q for (w, _), q in zip(omega, qs) if w)
    assert compare_quad_rational(M, A*B, D) >= 0
    terms = [(w*p, p) for (w, _), p in zip(tau, ps)]
    terms += [(w*q, q) for (w, _), q in zip(omega, qs)]
    N = lcm(*(w.denominator for w, _ in terms))
    # Z=exp(-N*(S_tau+S_omega)) is an exact rational, with integer powers.
    Z = F(1)
    for w, p in terms:
        e = N*w
        assert e.denominator == 1
        Z *= p**e.numerator
    lower_sign = compare_quad_rational(qpow(M, N, D), Z, D)
    assert lower_sign >= 0  # S_tau+S_omega >= -log M, source equation (11).
    mass = sum(w for w, _ in tau)*sum(w for w, _ in omega)
    assert Z >= mass**(-N)  # Source's upper entropy bound.
    rigid = d == 1 and all(w == 0 or vnorm2(a) == 1 for w, a in tau+omega)
    assert (lower_sign == 0) == rigid
    if rigid:
        assert sum(w for w, _ in tau) == sum(w for w, _ in omega) == D == Z == 1
    return lower_sign == 0


def main():
    assert mul(I, J) == K and mul(J, I) == scale(-1, K)
    fidelity = 0
    totals = []
    for d in (1, 2, 3):
        b = basis(d, False)
        pool = [tuple(ZERO for _ in range(d))] + b + basis(d, True)
        pool += [tuple(ONE for _ in range(d)), tuple(I if i % 2 else ONE for i in range(d))]
        for a, z, x in product(pool, repeat=3):
            assert buzano(a, z, x)
            fidelity += 1
        frames = [frame(d, t, r, k) for t in (False, True)
                  for r in (F(1), F(1, 2)) for k in (1, 2)]
        for fr in frames:
            parseval_matrix(fr, d)
        states = list(product((ONE, I, J, add(ONE, I)), repeat=d))
        good = eq = skipped = 0
        for tau, omega, x in product(frames, frames, states):
            ans = state_check(tau, omega, x, d)
            if ans is None:
                skipped += 1
            else:
                good += 1
                eq += ans
        totals.append((d, len(frames), len(states), good, eq, skipped))
    print('FIDELITY Theorem 2.1: exact quaternionic Buzano triples =', fidelity)
    for d, fs, xs, good, eq, skipped in totals:
        print(f'DIM={d} frames={fs} states={xs} admissible_pairs_states={good} '
              f'lower_equalities={eq} excluded_zero_coefficients={skipped}')
    # Null atoms are visible to admissibility but not to integrals.
    tau_zero = [(F(1), (ONE,)), (F(0), (ZERO,))]
    tau_small = [(F(1), (ONE,)), (F(0), (scale(F(1, 2), I),))]
    omega = [(F(1), (J,))]
    for tau in (tau_zero, tau_small):
        parseval_matrix(tau, 1)
    assert state_check(tau_zero, omega, (ONE,), 1) is None
    assert state_check(tau_small, omega, (ONE,), 1) is True
    print('NULL_ATOMS: zero frame vector destroys domain; nonzero subunit vector preserves equality')
    # Exact algebra used by the proposed general saturation route.
    factor_coefficients = [F(-1)*1, F(-1)*2+F(1)*1, F(1)*2]
    assert factor_coefficients == [-1, -1, 2]
    roots = (F(1), F(-1, 2))
    assert [r for r in roots if r >= 0] == [F(1)]
    for c in roots:
        assert 2*c*c-c-1 == 0
    print('SATURATION_IDENTITY: exact factor coefficients [-1,-1,2]; roots 1 and -1/2')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
