#!/usr/bin/env python3
"""Exact source-fidelity and singular-boundary residual checks.

All arithmetic is integer, rational or exact SymPy algebra. The finite checks
support the accompanying general mathematical argument; they do not certify an
infinite family or constitute Lean verification. No random or floating inputs.
"""
from itertools import product
import platform
import sympy as sp


def zero(matrix):
    assert all(sp.simplify(x) == 0 for x in matrix)


def eq(left, right):
    zero(left - right)


def source_binary():
    # Source: Conjecture 2's stated binary case, and Eqs. (3)-(4), i != j.
    s = sp.Symbol('s', positive=True)
    g = sp.Matrix([[1, s], [s, 1]])
    a = 1 - s  # valid for the source domain 0 < s < 1
    m = g - a * sp.eye(2)
    assert m.det() == 0
    eq(m / s, sp.ones(2))
    assert sp.simplify(sp.sqrt(s) * sp.sqrt(s) - s) == 0
    u, v = sp.symbols('u v', nonnegative=True)
    determinant = (sp.ones(2) - sp.diag(u, v)).det()
    assert sp.expand(determinant + u * (1 - v) + v) == 0
    # PSD and u,v in [0,1] imply u*(1-v)+v <= 0: v=0, then u=0.
    points = [sp.Rational(i, 8) for i in range(9)]
    feasible = []
    for x, y in product(points, repeat=2):
        r = sp.ones(2) - sp.diag(x, y)
        if r[0, 0] >= 0 and r[1, 1] >= 0 and r.det() >= 0:
            feasible.append((x, y))
    assert feasible == [(0, 0)]
    print('SOURCE_BINARY: symbolic residual=ones(2); det=-u*(1-v)-v; 81 exact grid points, only (0,0) feasible')


def family(n):
    b = sp.zeros(n)
    b[0, 0] = 1
    b[0, 1] = b[1, 1] = 1
    b[0, 2], b[1, 2] = 1, 2
    for j in range(3, n):
        b[0, j] = b[1, j] = b[j - 1, j] = 1
    weights = [1, 2, 5] + [3] * (n - 3)
    w = sp.diag(*weights)
    q = b.T * b
    assert list(q.diagonal()) == weights
    assert q.rank() == n - 1
    k = sp.Matrix([1, -2, 1] + [0] * (n - 3))
    zero(q * k)
    for i in range(n):
        for j in range(i):
            # Normalized overlap is Q_ij / sqrt(w_i w_j).
            assert 0 < q[i, j] ** 2 < weights[i] * weights[j]
    c = b[:2, :]
    bs = sp.symbols(f'b0:{n}', real=True)
    raw_residual = q - sp.diag(*(weights[i] * bs[i] for i in range(n)))
    charge = (k.T * raw_residual * k)[0]
    assert sp.expand(charge + bs[0] + 8 * bs[1] + 5 * bs[2]) == 0
    first_three_zero = {bs[i]: 0 for i in range(3)}
    box = raw_residual.subs(first_three_zero)
    eq(box, c.T * c + sp.diag(0, 0, 0, *(1 - 3 * bs[j] for j in range(3, n))))
    for j in range(3, n):
        x = sp.zeros(n, 1)
        x[j], x[1] = 1, -1
        assert sp.expand((x.T * box * x)[0] - (1 - 3 * bs[j])) == 0
        # An exact point within [0,1] but above the cap has a negative readout.
        assert (x.T * box * x)[0].subs(bs[j], sp.Rational(5, 12)) == -sp.Rational(1, 4)
    vertices = 0
    for bits in product((0, 1), repeat=n - 3):
        beta = sp.Matrix([0, 0, 0] + [sp.Rational(v, 3) for v in bits])
        r = q - w * sp.diag(*beta)
        factor = c.col_join(sp.diag(0, 0, 0, *(1 - v for v in bits)))
        eq(r, factor.T * factor)
        vertices += 1
    maximal_beta = sp.diag(0, 0, 0, *([sp.Rational(1, 3)] * (n - 3)))
    eq(q - w * maximal_beta, c.T * c)
    # G = p I + (1-p) H. Congruence by diag(sqrt(w_i)) is exact.
    # Reparameterize p=t/(1+t); G is congruent to (Q+tW)/(1+t).
    t = sp.Symbol('t', positive=True)
    r = q + t * w
    for size in range(1, n + 1):
        poly = sp.Poly(r[:size, :size].det(), t)
        assert all(a >= 0 for a in poly.all_coeffs())
        assert poly.LC() > 0 and poly.degree() == size
    for p in (sp.Rational(1, 5), sp.Rational(1, 2), sp.Rational(4, 5)):
        lifted = p * w + (1 - p) * q
        assert all(lifted[:size, :size].det() > 0 for size in range(1, n + 1))
        eq(lifted - p * w, (1 - p) * q)
        assert (lifted - p * w).det() == 0
        # Composition identity, after the same diagonal congruence.
        eq(lifted - w * (p * sp.eye(n) + (1 - p) * maximal_beta),
           (1 - p) * (q - w * maximal_beta))
        priors = [sp.Rational(i + 1, n * (n + 1) // 2) for i in range(n)]
        assert sum(priors) == 1 and all(x > 0 for x in priors)
        joint = p * sum(priors[3:]) / 3
        assert joint > 0
    print(f'FAMILY n={n}: rank={n-1}; pairwise 0<overlap<1; symbolic PD lift; exact residual box, {vertices} vertices; 3 rational p values')


def physical_instrument():
    # A complete physical n=4 instance, not just a PSD matrix witness.
    n, p = 4, sp.Rational(1, 2)
    b = sp.Matrix([[1, 1, 1, 1], [0, 1, 2, 1], [0, 0, 0, 1], [0, 0, 0, 0]])
    w = sp.diag(1, 2, 5, 3)
    sinv = sp.diag(1, 1 / sp.sqrt(2), 1 / sp.sqrt(5), 1 / sp.sqrt(3))
    phi = b * sinv
    h = phi.T * phi
    g = p * sp.eye(n) + (1 - p) * h
    # Cholesky is formed in rational coordinates before exact normalization.
    rational_g = p * w + (1 - p) * (b.T * b)
    psi = rational_g.cholesky().T * sinv
    eq(psi.T * psi, g)
    reciprocal = psi.inv()
    ks = []
    for i in range(n):
        ks.append(sp.sqrt(p) * phi[:, i] * reciprocal[i, :])
    k0 = sp.sqrt(1 - p) * phi * reciprocal
    eq(sum((k.T * k for k in ks), k0.T * k0), sp.eye(n))
    for j in range(n):
        eq(k0 * psi[:, j], sp.sqrt(1 - p) * phi[:, j])
        for i in range(n):
            eq(ks[i] * psi[:, j], sp.sqrt(p) * phi[:, j] if i == j else sp.zeros(n, 1))
    e4 = sp.diag(0, 0, 1, 0)
    e0 = sp.eye(n) - e4
    assert all(x >= 0 for x in e0.diagonal())
    for j in range(n):
        assert (phi[:, j].T * e4 * phi[:, j])[0] == (sp.Rational(1, 3) if j == 3 else 0)
    print('INSTRUMENT n=4,p=1/2: exact normalized inputs, Kraus completeness, zero error, outcome-independent pure outputs, residual POVM success=(0,0,0,1/3)')


def main():
    print(f'Python {platform.python_version()}; SymPy {sp.__version__}; exact arithmetic only')
    source_binary()
    for n in range(4, 9):
        family(n)
    physical_instrument()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
