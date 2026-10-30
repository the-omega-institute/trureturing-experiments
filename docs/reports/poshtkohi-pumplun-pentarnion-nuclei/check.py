#!/usr/bin/env python3
"""Exact checks for arXiv:2610.31485v1, Section 7 (nuclei conjunct).

Requires Python 3 and sympy 1.14.0. No floating-point computation is used.
The multiplication tensor is transcribed from Section 2, Table 2.
Finite tests support the proposed general proof; they are not its proof.
"""
from itertools import product
import sympy as s

a, b, c, d, gamma = s.symbols('a b c d gamma', real=True)
A, B, C, U, V, W = s.symbols('A B C U V W', real=True)
alpha, beta, zeta, delta = s.symbols('alpha beta zeta delta', real=True)
coeffs = (alpha, beta, zeta, delta)
zero = s.zeros(5, 1)
basis = [s.eye(5)[:, n] for n in range(5)]
one, i, j, k, ell = basis
tensor = {(p, q): zero.copy() for p, q in product(range(5), repeat=2)}
for p in range(5):
    tensor[0, p] = tensor[p, 0] = basis[p]
for p, square in enumerate((a, b, c, d), 1):
    tensor[p, p] = -square * one
for p, q, value in (
    (1, 2, A * (k + gamma * ell)),
    (1, 3, -B * (j + gamma * ell)),
    (2, 3, C * (i + gamma * ell)),
    (1, 4, -U * (j + k)),
    (2, 4, -V * (k + i)),
    (3, 4, -W * (i + j)),
):
    tensor[p, q], tensor[q, p] = value, -value


def mul(x, y):
    return sum((x[p] * y[q] * tensor[p, q]
                for p, q in product(range(5), repeat=2)
                if x[p] != 0 and y[q] != 0), zero.copy()).applyfunc(s.expand)


def assoc(x, y, z):
    return (mul(mul(x, y), z) - mul(x, mul(y, z))).applyfunc(s.expand)


def equal(x, y):
    if isinstance(x, s.MatrixBase):
        assert all(s.simplify(v) == 0 for v in x - y), (x, y)
    else:
        assert s.simplify(x - y) == 0, (x, y)


def conjugate(x):
    return s.diag(1, -1, -1, -1, -1) * x


x = s.Matrix([0, *coeffs])
E, F, G = a * V + b * U, a * W + c * U, b * W + c * V

# Fidelity: Theorem 3.9's normalized [x, ell, ell] and determinant 54.
normalized = {a: 1, b: 1, c: 1, d: 1,
              U: s.sqrt(2) / 2, V: s.sqrt(2) / 2, W: s.sqrt(2) / 2}
source_matrix = s.Matrix([[4, 1, 1], [1, 4, 1], [1, 1, 4]])
equal(assoc(x, ell, ell).subs(normalized),
      s.Matrix([0, *(source_matrix * s.Matrix(coeffs[:3]) / 2), 0]))
assert source_matrix.det() == 54
print('FIDELITY: Theorem 3.9 associator reproduced; determinant = 54')

# Independent symbolic checks, with the six radical coefficients initially free.
for u, v in product(basis, repeat=2):
    equal(conjugate(mul(u, v)), mul(conjugate(v), conjugate(u)))
    equal(assoc(one, u, v), zero)
    equal(assoc(u, one, v), zero)
    equal(assoc(u, v, one), zero)
equal(assoc(x, ell, i)[0], E * beta + F * zeta)
equal(assoc(x, ell, j)[0], E * alpha + G * zeta)
equal(assoc(x, ell, k)[0], F * alpha + G * beta)
equal(assoc(ell, i, i)[2], U * B)
middle_sum = assoc(i, x, j)[4] + assoc(j, x, k)[4] + assoc(k, x, i)[4]
equal(middle_sum, 2 * gamma * (C * U + B * V + A * W) * delta)
equal(assoc(i, x, j)[4].subs(delta, 0), gamma * A * (B * alpha - C * beta))
equal(assoc(j, x, k)[4].subs(delta, 0), gamma * C * (A * beta - B * zeta))
parameter = s.symbols('parameter', real=True)
equal(assoc(ell, x, i)[0].subs({alpha: C * parameter, beta: B * parameter,
                             zeta: A * parameter, delta: 0}),
      -a * (V * B + W * A) * parameter)
for u, q, v in zip((i, j, k), (a, b, c), coeffs):
    equal(assoc(u, x, ell)[4].subs(gamma, 0), -q * v)
equal(assoc(i, delta * ell, j)[0], E * delta)

# Squaring E=0 and clearing its nonzero denominators forces a=b.
equal(s.factor(a**2 * b * d * (b + c) - b**2 * a * d * (a + c)),
      a * b * c * d * (a - b))
# The F and G reductions are the corresponding variable permutations.
equal(s.factor(a**2 * c * d * (b + c) - c**2 * a * d * (a + b)),
      a * b * c * d * (a - c))
equal(s.factor(b**2 * c * d * (a + c) - c**2 * b * d * (a + b)),
      a * b * c * d * (b - c))
print('SYMBOLIC: scalar inclusion, conjugation, elimination and radical identities verified')

# Each row is an actual associator coordinate (or a sum of such coordinates).
# Nonzero 4-by-4 determinants force the four pure coordinates to vanish.
left_rows = [assoc(x, ell, u)[0] for u in (i, j, k)] + [assoc(x, i, i)[2]]
middle_rows_nonzero = [middle_sum, assoc(i, x, j)[4],
                      assoc(j, x, k)[4], assoc(ell, x, i)[0]]
middle_rows_zero = [assoc(u, x, ell)[4] for u in (i, j, k)] + [assoc(i, x, j)[0]]


def determinant(rows):
    return s.factor(s.Matrix([[s.expand(row).coeff(v) for v in coeffs]
                              for row in rows]).det())


left_det = determinant(left_rows)
middle_det_nonzero = determinant(middle_rows_nonzero)
middle_det_zero = s.factor(determinant(middle_rows_zero).subs(gamma, 0))
equal(left_det, 2 * E * F * G * U * B)
print('LEFT_MINOR =', left_det)
print('MIDDLE_MINOR_NONZERO_GAMMA =', middle_det_nonzero)
print('MIDDLE_MINOR_ZERO_GAMMA =', middle_det_zero)

tested = rejected = indefinite = gamma_zero = 0
values = (-2, -1, 1, 2)
gammas = (-2, -1, 0, 1, 2)
for aa, bb, cc, dd, gg in product(values, values, values, values, gammas):
    params = tuple(map(s.Integer, (aa, bb, cc, dd, gg)))
    aa, bb, cc, dd, gg = params
    denominators = (cc + gg**2 * dd, bb + gg**2 * dd, aa + gg**2 * dd,
                    bb + cc, aa + cc, aa + bb)
    if any(q == 0 for q in denominators):
        rejected += 1
        continue
    radicands = tuple(n / q for n, q in zip(
        (aa * bb, aa * cc, bb * cc, aa * dd, bb * dd, cc * dd), denominators))
    if any(q < 0 for q in radicands):
        rejected += 1
        continue
    assert all(q > 0 for q in radicands)
    substitution = dict(zip((a, b, c, d, gamma), params))
    substitution.update(zip((A, B, C, U, V, W), map(s.sqrt, radicands)))
    for det in (left_det, middle_det_zero if gg == 0 else middle_det_nonzero):
        value = s.simplify(det.subs(substitution))
        # Sympy's exact algebraic zero predicate; an unknown result fails.
        assert value.is_zero is False, (params, value)
    tested += 1
    indefinite += any(q < 0 for q in params[:4])
    gamma_zero += gg == 0
assert tested > 0 and indefinite > 0 and gamma_zero > 0
assert tested + rejected == 4**4 * 5
print(f'GRID: {tested} admissible tuples; {indefinite} indefinite; '
      f'{gamma_zero} with gamma=0; {rejected} outside source domain')
print('All three nuclei are the scalar line in every tested case: '
      'left and middle by nonzero minors, right by the checked conjugation identity.')
print('No Lean proof of the claim or its negation is included.')
print('ALL_CHECKS_OK')
