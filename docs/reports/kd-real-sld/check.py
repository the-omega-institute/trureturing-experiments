#!/usr/bin/env python3
"""Exact checks of the KD-real mixed-state SLD counterexample.

Requires Python 3 and SymPy. All inputs and assertions use exact algebra.
Source: Spriet et al., arXiv:2511.01996v2, Section 6.3 and Appendix E.
This is experimental evidence and symbolic identity checking, not a Lean proof.
"""
from itertools import product
import sympy as s


def normal(x):
    return s.factor(s.simplify(s.expand_complex(x)))


def matrix_normal(a):
    return a.applyfunc(normal)


def equal(a, b):
    if isinstance(a, s.MatrixBase):
        assert a.shape == b.shape
        assert matrix_normal(a - b) == s.zeros(*a.shape)
    else:
        assert normal(a - b) == 0


def comm(a, b):
    return a * b - b * a


def kd(a, b, c):
    """Left source KD distribution; columns of a and b are basis vectors."""
    u, v = a.H * b, a.H * c * b
    return s.Matrix(u.rows, u.cols, lambda x, k: normal(s.conjugate(u[x, k]) * v[x, k]))


def symbol(a, b, c):
    u, v = a.H * b, a.H * c * b
    assert all(normal(z) != 0 for z in u)
    return s.Matrix(u.rows, u.cols, lambda x, k: normal(v[x, k] / u[x, k]))


def real_matrix(a):
    assert all(normal(s.im(z)) == 0 for z in a)


def is_sld(rho, x, ell):
    equal(ell.H, ell)
    equal((rho * ell + ell * rho) / 2, -s.I * comm(x, rho))


def fidelity():
    """Reproduce the source's pure-state Proposition 6.5 and QFI normalization."""
    a = s.eye(2)
    b = s.Matrix([[1, 1], [1, -1]]) / s.sqrt(2)
    rho = s.diag(1, 0)
    x = s.Matrix([[0, 1], [1, 0]])
    ell = s.Matrix([[0, 2*s.I], [-2*s.I, 0]])
    equal(b.H * b, a)
    real_matrix(kd(a, b, rho))
    real_matrix(symbol(a, b, x))
    is_sld(rho, x, ell)
    equal(s.trace(rho * ell * ell), 4)
    equal(kd(a, b, ell), s.Matrix([[s.I, -s.I], [-s.I, s.I]]))
    pop = b.H * rho * b
    moment = b.H * ell * rho * b
    conditional = s.diag(*[normal(moment[k, k] / pop[k, k]) for k in range(2)])
    equal(conditional, s.diag(-2*s.I, 2*s.I))
    equal(conditional.H, -conditional)
    assert conditional != s.zeros(2)
    print('FIDELITY Proposition 6.5: pure qubit; QFI=4; L is not KD-real; conditional L is nonzero anti-Hermitian')


def graph_laws(p):
    """Exact modular phases for the monomial matrices, for every r,s,x in F_p^2."""
    g = list(product(range(p), repeat=2))

    def add(r, z):
        return tuple((r[i] + z[i]) % p for i in range(2))

    def neg(r):
        return tuple(-z % p for z in r)

    def j(r):
        return (-r[1] % p, r[0])

    def dot(r, z):
        return sum(r[i] * z[i] for i in range(2)) % p

    def act(r, x):
        y = add(x, r)
        return y, dot(j(r), y)  # W_r a_x = omega^phase a_y

    for r in g:
        assert dot(j(r), r) == 0
        for x in g:
            y, exponent = act(r, x)
            z, opposite = act(neg(r), y)
            assert z == x and (exponent + opposite) % p == 0
            # The KD phase of W_r is (Jr).x - k.r; opposite r conjugates it.
            for k in g:
                phase = (dot(j(r), x) - dot(k, r)) % p
                opposite_phase = (dot(j(neg(r)), x) - dot(k, neg(r))) % p
                assert (phase + opposite_phase) % p == 0
        for z in g:
            cocycle = -dot(j(z), r) % p
            for x in g:
                y, ez = act(z, x)
                y, er = act(r, y)
                target, ers = act(add(r, z), x)
                assert y == target and (ez + er - ers - cocycle) % p == 0
            if r != z:
                # W_r^* W_z has a nonzero translation, so every diagonal entry vanishes.
                shift = add(neg(r), z)
                assert all(add(x, shift) != x for x in g)
    shifts = {((u % p), (v % p)) for u in [-1, 1] for v in [-1, 1]}
    assert len(shifts) == 4
    # At displacement (1,1), [R,X] has (omega - omega^-1)/4.
    # A primitive odd-p root cannot satisfy omega=omega^-1, because p does not divide 2.
    assert 2 % p != 0
    print(f'GRAPH p={p}: all {p**6} monomial products, adjoints, opposite KD phases and trace orthogonality pass')


def symbolic_nine_dimensional_family():
    t = s.Symbol('t', real=True)
    p = 3
    g = list(product(range(p), repeat=2))
    index = {x: i for i, x in enumerate(g)}
    d = p*p
    omega = (-1 + s.sqrt(3)*s.I) / 2
    equal(omega**3, 1)
    assert normal(omega - 1) != 0
    a = s.eye(d)
    b = s.Matrix(d, d, lambda z, k: normal(omega**sum(g[z][i]*g[k][i] for i in range(2)) / p))
    equal(b.H * b, a)
    assert all(normal(s.conjugate(z)*z) == s.Rational(1, d) for z in b)

    def w(r):
        m = s.zeros(d)
        for x in g:
            y = tuple((x[i] + r[i]) % p for i in range(2))
            exponent = -r[1]*y[0] + r[0]*y[1]
            m[index[y], index[x]] = normal(omega**(exponent % p))
        return m

    pmat, qmat = w((1, 0)), w((0, 1))
    r, x = matrix_normal((pmat + pmat.H)/2), matrix_normal((qmat + qmat.H)/2)
    equal(pmat.H * pmat, a)
    equal(qmat.H * qmat, a)
    equal(pmat * qmat, omega**2 * qmat * pmat)
    equal(comm(r, x)[index[(1, 1)], index[(0, 0)]], (omega - s.conjugate(omega))/4)
    e = matrix_normal((a + pmat + pmat**2)/3)
    equal(e.H, e)
    equal(e*e, e)
    equal(s.trace(e), 3)
    equal(e*x*e, s.zeros(d))
    equal(e*x*x*e, e/2)

    rho = matrix_normal((a + t*r)/d)
    alpha, beta = (1+t)/9, (1-t/2)/9
    equal(rho, alpha*e + beta*(a-e))
    equal(s.trace(rho), 1)
    # Exact positivity certificates on 0<t<1: alpha>1/9 and beta>1/18.
    equal(alpha-s.Rational(1, 9), t/9)
    equal(beta-s.Rational(1, 18), (1-t)/18)
    equal(rho.H, rho)
    equal(s.trace(rho*rho), (1+t*t/2)/9)
    assert matrix_normal(comm(x, rho)) != s.zeros(d)
    qrho = kd(a, b, rho)
    equal(qrho, s.Matrix(d, d, lambda z, k: normal((1 + t*(omega**((g[z][1]-g[k][0]) % p) + omega**(-(g[z][1]-g[k][0]) % p))/2) / d**2)))
    kd_values = {normal(v) for v in qrho}
    assert kd_values == {normal((1+t)/81), normal((2-t)/162)}
    # These are >1/81 and >1/162 respectively on the stated open interval.
    real_matrix(symbol(a, b, rho))
    real_matrix(symbol(a, b, x))
    ell = matrix_normal(-6*s.I*t/(4+t)*comm(x, e))
    is_sld(rho, x, ell)
    real_matrix(kd(a, b, ell))
    real_matrix(symbol(a, b, ell))
    qfi = normal(s.trace(rho*ell*ell))
    equal(qfi, 3*t*t/(4+t))
    # Positive QFI for 0<t<1: its numerator is positive and denominator exceeds 4.
    for basis in [a, b]:
        pop, deriv = matrix_normal(basis.H*rho*basis), matrix_normal(basis.H*(-s.I*comm(x, rho))*basis)
        for k in range(d):
            equal(pop[k, k], s.Rational(1, d))
            equal(deriv[k, k], 0)
    print('FAMILY p=3, all real 0<t<1: exact symbolic density, strict positivity certificates, KD positivity, SLD and QFI=3*t^2/(4+t)')
    print('FIDELITY Proposition 6.4: computational/Fourier measurement marginals=1/9 and tangent diagonals=0; both classical FI vanish')
    rh, lh = rho.subs(t, s.Rational(1, 2)), ell.subs(t, s.Rational(1, 2))
    equal(rh, (a+e)/12)
    equal(lh, -2*s.I*comm(x, e)/3)
    equal(qfi.subs(t, s.Rational(1, 2)), s.Rational(1, 6))
    assert sorted({normal(v) for v in qrho.subs(t, s.Rational(1, 2))}) == [s.Rational(1, 108), s.Rational(1, 54)]
    assert all(not z.has(s.Float) for m in [rho, ell, qrho, rh, lh] for z in m)
    print('MEMBER p=3,t=1/2: rho=(I+E)/12, eigenvalues 1/6 (x3), 1/12 (x6); KD entries {1/54,1/108}; QFI=1/6')


def main():
    print('Exact arithmetic: Python integers and SymPy '+s.__version__+'; no numerical approximations')
    fidelity()
    for p in [3, 5, 7]:
        graph_laws(p)
    symbolic_nine_dimensional_family()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
