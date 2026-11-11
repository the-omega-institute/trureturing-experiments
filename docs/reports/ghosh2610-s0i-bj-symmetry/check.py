#!/usr/bin/env python3
"""Exact checks supporting Ghosh's S0I/Birkhoff--James symmetry preregistration.

Only integers and fractions are used. The checks do not prove the universal
conjecture. The analytic norm bound uses Holder's inequality, as explained in
README.md; the program verifies the exact algebra and its explicit witnesses.
"""

from dataclasses import dataclass
from fractions import Fraction as F
from itertools import product


def eye(n):
    return tuple(tuple(F(i == j) for j in range(n)) for i in range(n))


def add(a, b):
    return tuple(tuple(x + y for x, y in zip(r, s)) for r, s in zip(a, b))


def scale(c, a):
    return tuple(tuple(c * x for x in r) for r in a)


def transpose(a):
    return tuple(zip(*a))


def mul(a, b):
    return tuple(tuple(sum(x * y for x, y in zip(r, c)) for c in zip(*b)) for r in a)


def apply(a, v):
    return tuple(sum(x * y for x, y in zip(r, v)) for r in a)


def dot(u, v):
    return sum(x * y for x, y in zip(u, v))


def norm_power(v, p):
    return sum(abs(x) ** p for x in v)


def dual_gradient(v, p):
    return tuple((1 if x > 0 else -1 if x < 0 else 0) * abs(x) ** (p - 1) for x in v)


def fidelity_example_46():
    t = ((F(0), F(0), F(1)), (F(0), F(2), F(0)), (F(1), F(0), F(0)))
    inv = ((F(0), F(0), F(1)), (F(0), F(1, 2), F(0)), (F(1), F(0), F(0)))
    assert mul(t, inv) == mul(inv, t) == eye(3)
    e1, e2 = (F(1), F(0), F(0)), (F(0), F(1), F(0))
    assert dot(e1, apply(t, e1)) == dot(e1, apply(inv, e1)) == 0
    gram = mul(transpose(t), t)
    assert gram == ((F(1), F(0), F(0)), (F(0), F(4), F(0)), (F(0), F(0), F(1)))
    # 4||v||_2^2 - ||Tv||_2^2 = 3(v_1^2 + v_3^2): M_T = {+/- e_2}.
    assert add(scale(4, eye(3)), scale(-1, gram)) == (
        (F(3), F(0), F(0)), (F(0), F(0), F(0)), (F(0), F(0), F(3)))
    shifted = add(t, scale(F(-1, 2), eye(3)))
    row = ((F(1), F(0), F(1)),)
    # (9/4)I - (T-I/2)^*(T-I/2) = row^* row, an exact PSD certificate.
    assert add(scale(F(9, 4), eye(3)), scale(-1, mul(transpose(shifted), shifted))) == mul(transpose(row), row)
    assert norm_power(apply(shifted, e2), 2) == F(9, 4) < 4
    print('FIDELITY Example 4.6: inverse, 0 in both numerical ranges, M_T={+/-e2}, ||T-I/2||=3/2<2 verified')


@dataclass(frozen=True)
class Gaussian:
    re: F = F(0)
    im: F = F(0)

    def __add__(self, other):
        return Gaussian(self.re + other.re, self.im + other.im)

    def __neg__(self):
        return Gaussian(-self.re, -self.im)

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        return Gaussian(self.re * other.re - self.im * other.im,
                        self.re * other.im + self.im * other.re)

    def __truediv__(self, other):
        d = other.re ** 2 + other.im ** 2
        assert d != 0
        return self * Gaussian(other.re / d, -other.im / d)


ZERO, ONE = Gaussian(), Gaussian(F(1))


def idempotent_product(a, b):
    """Multiply a_P P + a_I I modulo P^2=P, over exact Gaussian rationals."""
    ap, ai = a
    bp, bi = b
    return (ap * bp + ap * bi + ai * bp, ai * bi)


def shifted_idempotents():
    scalars = sorted({F(n, d) for n in range(-4, 5) for d in range(1, 4)})
    zs = [Gaussian(a, b) for a, b in product(scalars, repeat=2)]
    n = 0
    for z in zs:
        if z in (ZERO, ONE):
            continue
        s = (ONE, -z)
        inv = (ONE / (z * (ONE - z)), -ONE / z)
        assert idempotent_product(s, inv) == (ZERO, ONE)
        assert idempotent_product(inv, s) == (ZERO, ONE)
        # Evaluation at a supporting pair: inverse value zero forces g(Pv)=1-z.
        assert inv[0] * (ONE - z) + inv[1] == ZERO
        n += 1
    resolvent_count = 0
    for t in scalars:
        if t <= 0:
            continue
        tg = Gaussian(t)
        r = (-tg / (ONE + tg), ONE)
        expansion = (tg, ONE)
        assert idempotent_product(r, expansion) == (ZERO, ONE)
        assert idempotent_product(expansion, r) == (ZERO, ONE)
        assert r[0] + ONE == ONE / (ONE + tg)
        resolvent_count += 1
    print(f'IDENTITIES: both shifted-projection inverse identities at {n} real/complex rational shifts and {resolvent_count} positive rational resolvents')


def lp_contrapositive():
    total = 0
    for p in range(3, 9):
        u = (F(1), F(2))
        y = (F(2 ** (p - 1)), F(-1))
        ju = dual_gradient(u, p)
        c = norm_power(u, p)
        assert dot(ju, y) == 0  # Holder gives u perpendicular_B y over R and C.
        assert dot(dual_gradient(y, p), u) != 0
        alpha = next(F(-1, 2 ** k) for k in range(1, 20)
                     if norm_power(tuple(a + F(-1, 2 ** k) * b for a, b in zip(y, u)), p) < norm_power(y, p))
        projector = tuple(tuple(a * b / c for b in ju) for a in u)
        assert mul(projector, projector) == projector
        # j(u)_i = sign(u_i)|u_i|^(p-1) certifies the Holder equality exponents.
        assert all(abs(a) ** p == abs(b) ** (p * (p - 1)) for a, b in zip(ju, u))
        v = (F(2), F(-1))
        jv = dual_gradient(v, p)
        z = dot(jv, apply(projector, v)) / norm_power(v, p)
        assert z == -F((2 ** (p - 1) - 2) ** 2, (1 + 2 ** p) ** 2) < 0
        t = add(projector, scale(-z, eye(2)))
        inv = add(scale(1 / (z * (1 - z)), projector), scale(-1 / z, eye(2)))
        assert mul(t, inv) == mul(inv, t) == eye(2)
        assert dot(jv, apply(t, v)) == 0
        # ||P||=1 implies Re g(Pw)<=1 for every supporting pair (w,g).
        # Thus Re g(T^-1 w)>=(1-(1-z))/(z(1-z))=1/(1-z)>0.
        lower = 1 / (1 - z)
        assert (1 - (1 - z)) / (z * (1 - z)) == lower > 0
        assert 1 - z > 1
        for a, b in product(range(-3, 4), repeat=2):
            if a == b == 0:
                continue
            w = (F(a), F(b))
            value = dot(dual_gradient(w, p), apply(inv, w)) / norm_power(w, p)
            assert value >= lower
            total += 1
        print(f'LP p={p}: u={u}, y={y}, reverse fails at alpha={alpha}, z={z}, inverse real-part lower bound={lower}')
        if p == 4:
            assert t == ((F(53, 289), F(136, 289)), (F(34, 289), F(308, 289)))
            print('EXPLICIT p=4: T=(1/289)[[53,136],[34,308]], z=-36/289, lower=289/325; valid over R and C')
    print(f'CONTRAPOSITIVE: 6 exact lp families and {total} nonzero integer test vectors verified; universal exclusion uses the stated analytic norm bound')


def main():
    fidelity_example_46()
    shifted_idempotents()
    lp_contrapositive()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
