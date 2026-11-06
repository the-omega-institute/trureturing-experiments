#!/usr/bin/env python3
"""Exact source-fidelity and collision checks for BCMN Conjecture 13.

No floating-point arithmetic or numerical eigensolver is used. Source:
arXiv:1302.4618v3, Section 2.2; Candan--Kutay--Ozaktas (2000),
equations (16), (22), and Table IV. The Candan sign is negative;
the displayed three-dimensional BCMN example uses its complex conjugate.
"""
import sys
import sympy as s
from sympy.polys.matrices import DomainMatrix

I = s.I
R = s.Rational
K = s.QQ.algebraic_field(s.sqrt(2), s.sqrt(3))


def exact(z):
    """Canonical arithmetic in Q(sqrt(2),sqrt(3),i)."""
    z = s.expand(z, complex=True)
    re, im = z.as_real_imag()
    return K.to_sympy(K.from_sympy(re)) + I * K.to_sympy(K.from_sympy(im))


def mat_exact(a):
    return a.applyfunc(exact)


def equal(a, b):
    d = a - b
    if isinstance(d, s.MatrixBase):
        assert all(exact(z) == 0 for z in d), d
    else:
        assert exact(d) == 0, d


def norm_sq(z):
    return exact(s.conjugate(z) * z)


def hermitian_basis(n):
    result = []
    for j in range(n):
        e = s.zeros(n)
        e[j, j] = 1
        result.append(e)
    for j in range(n):
        for k in range(j + 1, n):
            e = s.zeros(n)
            e[j, k] = e[k, j] = 1
            result.append(e)
            e = s.zeros(n)
            e[j, k], e[k, j] = I, -I
            result.append(e)
    return result


def labels(m):
    return list(range(m)) if m % 2 else list(range(m - 1)) + [m]


def source_s(m):
    return s.Matrix(m, m, lambda j, k:
        (2 * s.cos(2 * s.pi * R(j, m)) - 4 if j == k else 0)
        + int((j + 1) % m == k) + int((k + 1) % m == j))


def fidelity_m3():
    """Audit the displayed M=3 transform and expose its false HMW assertion."""
    r = s.sqrt(3)
    p0 = s.Matrix([[3+r, r, r],
                   [r, (3-r)/2, (3-r)/2],
                   [r, (3-r)/2, (3-r)/2]]) / 6
    p2 = s.Matrix([[3-r, -r, -r],
                   [-r, (3+r)/2, (3+r)/2],
                   [-r, (3+r)/2, (3+r)/2]]) / 6
    p1 = s.Matrix([[0, 0, 0], [0, 1, -1], [0, -1, 1]]) / 2
    projectors = [p0, p1, p2]
    equal(sum(projectors, s.zeros(3)), s.eye(3))
    for j, p in enumerate(projectors):
        equal(p*p, p)
        equal(p.T, p)
        equal(s.trace(p), 1)
        for k, q in enumerate(projectors):
            if j != k:
                equal(p*q, s.zeros(3))
    # The four literal phases in the displayed BCMN formula, with its + sign.
    eighth = (1 + I) / s.sqrt(2)
    phases = [(1, 1), (I, eighth), (-1, I), (-I, (-1 + I)/s.sqrt(2))]
    fs = [mat_exact(p0 + z2*p2 + z1*p1) for z2, z1 in phases]
    dft = s.Matrix(3, 3, lambda j, k:
                   s.expand_complex(s.exp(2*s.pi*I*R(j*k, 3))) / s.sqrt(3))
    equal(fs[2], dft)
    basis = hermitian_basis(3)
    rows = []
    for f in fs:
        equal(f.conjugate().T*f, s.eye(3))
        for j in range(3):
            v = f[:, j]
            rows.append([exact((v.conjugate().T*b*v)[0]) for b in basis])
    a = s.Matrix(rows)
    dm = DomainMatrix.from_Matrix(a).convert_to(K)
    rank = dm.rank()
    null = dm.nullspace().to_Matrix()
    assert rank == 6 and null.rows == 3
    h = mat_exact(sum((null[0, j]*b for j, b in enumerate(basis)), s.zeros(3)))
    # Normalize merely to make the certificate and printed determinant compact.
    first = next(z for z in h if z != 0)
    h = mat_exact(h / first)
    for f in fs:
        equal(s.diag(*(exact((f[:, j].conjugate().T*h*f[:, j])[0])
                       for j in range(3))), s.zeros(3))
    det = exact(h.det())
    assert h.rank() == 2 and det == 0
    plus = s.Matrix([1/s.sqrt(2), R(1, 2), -R(1, 2)])
    minus = s.Matrix([1/s.sqrt(2), -R(1, 2), R(1, 2)])
    equal((plus.conjugate().T*plus)[0], 1)
    equal((minus.conjugate().T*minus)[0], 1)
    equal((plus.conjugate().T*minus)[0], 0)
    expected = [[R(1, 2), R(1, 4), R(1, 4)], [R(1, 3)]*3,
                [R(1, 6), R(5, 12), R(5, 12)], [R(1, 3)]*3]
    for f, probabilities in zip(fs, expected):
        assert [norm_sq(z) for z in f.conjugate().T*plus] == probabilities
        assert [norm_sq(z) for z in f.conjugate().T*minus] == probabilities
    print('SOURCE DISCREPANCY M=3: four displayed bases are unitary; analysis rank=6, not 8')
    print('SOURCE DISCREPANCY M=3: rank-two nullspace generator =', h.tolist())
    print('SOURCE DISCREPANCY M=3: orthogonal unit collision gives', expected)
    print('SOURCE DISCREPANCY M=3: claimed successful HMW test does not reproduce')


def source_fidelity():
    """Reproduce the actual successful 3x8 HMW example in source equation (8)."""
    phi = s.Matrix([[2, 1, 1, 0, 0, 0, 1, I],
                    [-1, 0, 0, 1, 1, -1, -2, 2],
                    [0, 1, -1, 1, -1, 2*I, I, -1]])
    basis = hermitian_basis(3)
    rows = [[exact((phi[:, j].conjugate().T*b*phi[:, j])[0])
             for b in basis] for j in range(phi.cols)]
    a = s.Matrix(rows)
    assert a.rank() == 8
    kernel = a.nullspace()
    assert len(kernel) == 1
    h = sum((z*b for z, b in zip(kernel[0], basis)), s.zeros(3))
    h = mat_exact(h / next(z for z in h if z != 0))
    det = exact(h.det())
    assert det != 0
    for j in range(phi.cols):
        equal((phi[:, j].conjugate().T*h*phi[:, j])[0], 0)
    print('FIDELITY source 3x8 example: analysis rank=8, nullity=1')
    print('FIDELITY source 3x8 example: normalized kernel generator=', h.tolist())
    print('FIDELITY source 3x8 example: determinant=', det, '; HMW returns INJECTIVE')


def witness_m4():
    r = s.sqrt(2)
    u0 = s.Matrix([2+r, r, 2-r, r]) / 4
    u1 = s.Matrix([0, 1, 0, -1]) / r
    u2 = s.Matrix([1, -1, -1, -1]) / 2
    u4 = s.Matrix([2-r, -r, 2+r, -r]) / 4
    u = s.Matrix.hstack(u0, u1, u2, u4)
    eigs = [2*r - 4, -4, -4, -2*r - 4]
    equal(u.T*u, s.eye(4))
    smat = source_s(4)
    for j, k in enumerate(labels(4)):
        equal(smat*u[:, j], eigs[j]*u[:, j])
        for t in range(4):
            equal(u[(-t) % 4, j], (-1)**k*u[t, j])
    assert eigs[0] > eigs[2] > eigs[3]  # even parity ordering, odd block is 1D
    plus, minus = (u0 + I*u4)/r, (u0 - I*u4)/r
    equal((plus.conjugate().T*plus)[0], 1)
    equal((minus.conjugate().T*minus)[0], 1)
    equal((plus.conjugate().T*minus)[0], 0)
    assert all(norm_sq(z) != 0 for z in plus)
    expected = [R(3, 8), R(1, 8), R(3, 8), R(1, 8)]
    h = mat_exact(plus*plus.conjugate().T - minus*minus.conjugate().T)
    equal(h, I*(u4*u0.T-u0*u4.T))
    equal(h*h*h, h)
    equal(s.trace(h), 0)
    equal(s.trace(h*h), 2)
    assert h.rank() == 2
    for a in [R(0), R(1, 2), R(1), R(3, 2)]:
        diagonal = s.diag(*(s.expand_complex(s.exp(-I*s.pi*k*a/2)) for k in labels(4)))
        f = mat_exact(u*diagonal*u.T)
        equal(f.conjugate().T*f, s.eye(4))
        for state in [plus, minus]:
            actual = [norm_sq(z) for z in f.conjugate().T*state]
            assert actual == expected, (a, actual)
        equal(s.diag(*(exact((f[:, j].conjugate().T*h*f[:, j])[0])
                       for j in range(4))), s.zeros(4))
    dft = s.Matrix(4, 4, lambda j, k: (-I)**(j*k)/2)
    equal(u*s.diag(*(s.exp(-I*s.pi*k/2) for k in labels(4)))*u.T, dft)
    print('WITNESS M=4: labels=[0,1,2,4]; exact S eigenvalues=', eigs)
    print('WITNESS M=4: orthogonal unit states, all coordinates nonzero')
    print('WITNESS M=4: all four column-basis intensity vectors=', expected)
    print('WITNESS M=4: Hermitian certificate rank=2, eigenvalues +1 and -1')


def algebraic_family():
    # This identity holds for every real coordinate pair, not just sampled vectors.
    t, v = s.symbols('t v', real=True)
    for sign in [s.Integer(-1), s.Integer(1)]:
        z1, z2 = (t+I*sign*v)/s.sqrt(2), (t-I*sign*v)/s.sqrt(2)
        assert s.expand(s.conjugate(z1)*z1 - s.conjugate(z2)*z2) == 0
        assert s.expand(s.conjugate(z1)*z1) == (t*t+v*v)/2
    # Inner products in the real orthonormal two-mode coordinates.
    cplus = s.Matrix([1, I])/s.sqrt(2)
    cminus = s.Matrix([1, -I])/s.sqrt(2)
    equal((cplus.conjugate().T*cminus)[0], 0)
    equal((cplus.conjugate().T*cplus)[0], 1)
    equal((cminus.conjugate().T*cminus)[0], 1)
    p = s.symbols('p', integer=True)
    q = s.symbols('q', positive=True, integer=True)
    a0 = s.symbols('a0', real=True)
    phase = s.simplify(s.exp(I*s.pi*(2*q)*(a0+p/q)/2)
                       * s.exp(-I*s.pi*(2*q)*a0/2))
    assert s.simplify(phase - (-1)**p) == 0
    assert s.simplify(phase**2 - 1) == 0
    # These finite label checks supplement the general membership argument;
    # they do not certify arbitrary-dimensional source eigenbasis existence.
    assert all(0 in labels(m) and 4 in labels(m) for m in range(4, 129))
    cases = 0
    for denominator in range(1, 33):
        lower = max(4, 2*denominator)
        for m in range(lower, lower+10):
            assert 0 in labels(m) and 2*denominator in labels(m)
            cases += 1
    print('FAMILY: arbitrary real coordinates give equal signed imaginary-coherence intensities')
    print('FAMILY: symbolic affine-grid relative phase = (-1)^p, for all integer p and Q>=1')
    print('FAMILY: label checks M=4..128; affine-grid boundary/nearby cases=', cases)
    print('BOUNDARY: M=3 has no mode 4; general eigenbasis existence is sourced, not numerically inferred')


if __name__ == '__main__':
    print('Python', sys.version.split()[0], '; SymPy', s.__version__, '; exact arithmetic only')
    source_fidelity()
    fidelity_m3()
    witness_m4()
    algebraic_family()
    print('ALL_CHECKS_OK')
