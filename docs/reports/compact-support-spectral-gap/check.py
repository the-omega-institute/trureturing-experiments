#!/usr/bin/env python3
"""Exact helper checks for Kerner--Taufer Remark 2.10; no numerical spectral proof.

The continuum theorem is a general analytic target.  This program checks source
identities, formal quantization coefficients, and every member of a specified
finite family of genuine positive piecewise-constant barriers.  It does not
certify the infinite-dimensional theorem or its asserted open status.
"""
from fractions import Fraction as F
from itertools import product
import hashlib
from pathlib import Path
import sympy as S


def matrix(a, b, c, d):
    return ((F(a), F(b)), (F(c), F(d)))


def mul(t, u):
    return tuple(tuple(sum(t[i][k] * u[k][j] for k in range(2))
                       for j in range(2)) for i in range(2))


def determinant(t):
    return t[0][0] * t[1][1] - t[0][1] * t[1][0]


def free(a):
    return matrix(1, a, 0, 1)


def barrier(kappa, q):
    # Actual v=kappa**2 on length log(q)/kappa, q=exp(kappa*length).
    ch, sh = (F(q) + 1 / F(q)) / 2, (F(q) - 1 / F(q)) / 2
    return matrix(ch, sh / kappa, kappa * sh, ch)


def inverse_response(t):
    a, b = t[0]
    c, d = t[1]
    assert determinant(t) == 1
    return matrix(d / c, 1 / c, 1 / c, a / c)


def check_transfer(t):
    a, b = t[0]
    c, d = t[1]
    assert determinant(t) == 1 and a >= 1 and d >= 1 and b > 0 and c > 0
    m = matrix(a / b, -1 / b, -1 / b, d / b)
    n = inverse_response(t)
    assert mul(m, n) == matrix(1, 0, 0, 1)
    assert determinant(m) == c / b > 0 and m[0][0] > 0
    # Squared separation of the two eigenvalues of M^{-1}, with no radicals.
    separation_squared = (n[0][0] + n[1][1]) ** 2 - 4 * determinant(n)
    assert separation_squared == ((a - d) ** 2 + 4) / c ** 2 > 0
    # The coefficient C_v/pi^2 is 16*sqrt(separation_squared).
    for padding in (F(0), F(1, 3), F(2)):
        nt = inverse_response(mul(free(padding), mul(t, free(padding))))
        assert nt == matrix(n[0][0] + padding, n[0][1],
                            n[1][0], n[1][1] + padding)
    return separation_squared


def fidelity():
    x, y, length = S.symbols('x y L', real=True, positive=True)
    phi = S.Function('phi')
    potential = S.Function('v')
    transformed = S.sqrt(length) * phi(length * y)
    gl = -S.diff(transformed, y, 2) + length ** 2 * potential(length * y) * transformed
    uhl = S.sqrt(length) * (-S.diff(phi(x), x, 2) + potential(x) * phi(x)).subs(x, length * y)
    assert S.simplify(gl - length ** 2 * uhl) == 0
    # Proposition 2.4's scaling: epsilon_j=L^{-2} lambda_j.
    for n in (1, 2):
        u = S.sin(n * S.pi * (x + length / 2) / length)
        energy = (n * S.pi / length) ** 2
        assert S.simplify(-S.diff(u, x, 2) - energy * u) == 0
        assert S.simplify(u.subs(x, -length / 2)) == 0
        assert S.simplify(u.subs(x, length / 2)) == 0
    assert S.simplify((2 * S.pi / length) ** 2 - (S.pi / length) ** 2
                      - 3 * S.pi ** 2 / length ** 2) == 0
    print('FIDELITY: source unitary-scaling identity and free gap 3*pi**2/L**2: OK')


def formal_quantization():
    x = S.Symbol('x')  # x=1/ell, formal expansion variable.
    b = S.Symbol('b', positive=True)
    d = S.Symbol('d', real=True)  # Arbitrary first derivative of b(E) at zero.
    k = S.pi * x - S.pi / b * x ** 2 + S.pi / b ** 2 * x ** 3
    residual = S.series(k / x + S.atan(k / (b + d * k ** 2)) - S.pi,
                        x, 0, 3).removeO()
    assert S.simplify(residual) == 0
    e = S.series(k ** 2, x, 0, 5).removeO().expand()
    assert e == S.pi ** 2 * x ** 2 - 2 * S.pi ** 2 / b * x ** 3 + 3 * S.pi ** 2 / b ** 2 * x ** 4
    bm, bp = S.symbols('b_minus b_plus', positive=True)
    splitting = S.expand(e.subs(b, bp) - e.subs(b, bm))
    assert S.simplify(splitting.coeff(x, 3) - 2 * S.pi ** 2 * (1 / bm - 1 / bp)) == 0
    z, r = S.symbols('z R')
    # ell=L/2-R and z=1/L: the cubic conversion contributes factor 8.
    assert S.series((z / (S.Rational(1, 2) - r * z)) ** 3, z, 0, 4).removeO() == 8 * z ** 3
    print('QUANTIZATION: symbolic root coefficients and factor 16*pi**2: OK')


def small_barriers():
    types = tuple(product((1, 2, 3), (2, 3)))
    count = 0
    asymmetric = 0
    for segments in (1, 2, 3):
        for blocks in product(types, repeat=segments):
            for gaps in product((F(0), F(1, 2)), repeat=segments - 1):
                t = matrix(1, 0, 0, 1)
                for j, (kappa, q) in enumerate(blocks):
                    if j:
                        t = mul(free(gaps[j - 1]), t)
                    t = mul(barrier(kappa, q), t)
                check_transfer(t)
                count += 1
                asymmetric += int(t[0][0] != t[1][1])
    assert count == 942 and asymmetric > 0
    print(f'BARRIERS: all {count} cases (1-3 blocks, kappa=1,2,3, exp(kappa*width)=2,3, gaps=0,1/2); {asymmetric} asymmetric: OK')
    # Source's symmetric step special case: C_v=32*pi^2/(kappa*sinh(kappa*width)).
    for kappa, q in types:
        t = barrier(kappa, q)
        a, b = t[0]
        c, _ = t[1]
        bm, bp = (a - 1) / b, (a + 1) / b
        assert 0 < bm < bp
        assert 16 * (1 / bm - 1 / bp) == 32 / c
    print('STEP: all six exact symmetric-step coefficients agree with 32*pi**2/(kappa*sinh(kappa*width)): OK')


def main():
    fidelity()
    formal_quantization()
    small_barriers()
    print('SCOPE: exact helper identities and finite barrier family; not a proof of the universal analytic claim')
    print('CHECK_SHA256: ' + hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
