#!/usr/bin/env python3
"""Exact checks for Shapiro, arXiv:2605.14817v2, Conjecture 8.1.

All polynomial computations take place over ZZ (or a polynomial ring over ZZ).
The dimension readings are integer arithmetic, not a computation of the
irreducible components. Corollary 6.6 supplies the general spin identity;
the family argument and geometric justification are stated in README.md.
"""
import sympy as sp

L, t, w = sp.symbols('lambda t w')
u, v, z, m = sp.symbols('u v z m')


def continuant(a, c):
    assert len(c) + 1 == len(a)
    previous, current = sp.Poly(1, L, t), sp.Poly(L + a[0], L, t)
    for k in range(1, len(a)):
        previous, current = current, (
            sp.Poly(L + a[k], L, t) * current
            - sp.Poly(t * c[k - 1], L, t) * previous)
    return current


def spin_data(n):
    return ([n + 1 - 2*i for i in range(1, n + 1)],
            [i*(n-i) for i in range(1, n)])


def spin_product(n):
    qs = range(1, n, 2) if n % 2 == 0 else range(2, n, 2)
    return sp.Poly((L if n % 2 else 1) *
                   sp.prod(L**2-q*q*(1+t) for q in qs), L, t)


def equal(p, q):
    assert (sp.Poly(p, L, t) - sp.Poly(q, L, t)).is_zero


def main():
    print('Arithmetic: exact integer polynomial arithmetic; no numerical approximation')
    print('Source: Boris Shapiro, arXiv:2605.14817v2, Corollary 6.6 and Conjecture 8.1')
    # Source fidelity: reproduce both parity formulas, not merely the n=16 witness.
    for n in range(2, 25):
        a, c = spin_data(n)
        p = continuant(a, c)
        equal(p, spin_product(n))
        assert len(set(a)) == n and min(c) > 0
        if n in (4, 6):
            print(f'FIDELITY n={n}: P_n = {sp.factor(p.as_expr())}')
            # A separate determinant algorithm cross-checks the recurrence. The
            # asymmetric pencil has off-diagonal products w^2 c_i, so has the
            # same continuant as the symmetric source pencil with b_i^2=c_i.
            matrix = sp.diag(*(L + ai for ai in a))
            for i, ci in enumerate(c):
                matrix[i, i+1], matrix[i+1, i] = w*ci, w
            det = sp.Poly(sp.expand(matrix.det(method='domain-ge')), L, w)
            expected = sp.Poly(p.as_expr().subs(t, w*w), L, w)
            assert det == expected
            print(f'INDEPENDENT_DETERMINANT n={n}: exact equality')
    print('FIDELITY Corollary 6.6: all n=2,...,24 passed (both parities)')
    # Fixed counterexample: grouping several factors is allowed by Section 7.1.
    n = 16
    a, c = spin_data(n)
    indices = list(range(5, 13))  # source's one-based indices
    f = sp.Poly(sp.prod(L**2-q*q*(1+t) for q in (1, 3, 5, 7)), L, t)
    g = sp.Poly(sp.prod(L**2-q*q*(1+t) for q in (9, 11, 13, 15)), L, t)
    p = continuant(a, c)
    assert p == f*g
    assert sp.Poly(f.as_expr().subs(t, 0), L) == sp.Poly(
        sp.prod(L+a[i-1] for i in indices), L)
    assert f.degree(L) == 8 and g.degree(L) == 8
    print(f'WITNESS n=16: a={a}')
    print(f'WITNESS n=16: c={c}')
    print(f'WITNESS S={indices}: P_16=FG, F(lambda,0)=product_(i in S)(lambda+a_i)')
    print(f'WITNESS monomial counts: P={len(p.terms())}, F={len(f.terms())}, G={len(g.terms())}')
    kappa = (8*(16-8)+1)//2
    ambient = 2*16-1
    assert (ambient, kappa, kappa-ambient) == (31, 32, 1)
    print(f'DIMENSION n=16,d=8: ambient={ambient}, predicted_codimension={kappa}, excess={kappa-ambient}')
    # An all-m symbolic arithmetic identity, together with m>=4, proves the
    # numerical obstruction once Corollary 6.6 has supplied the grouped factor.
    difference = 2*m*m-(8*m-1)
    assert sp.expand(difference-(2*m*(m-4)+1)) == 0
    print('FAMILY exact identity: kappa(4m,2m)-(8m-1)=2m(m-4)+1 >= 1 for m>=4')
    for mm in range(4, 9):
        n = 4*mm
        a, c = spin_data(n)
        selected = range(mm+1, 3*mm+1)
        f = sp.Poly(sp.prod(L**2-q*q*(1+t) for q in range(1, 2*mm, 2)), L, t)
        g = sp.Poly(sp.prod(L**2-q*q*(1+t) for q in range(2*mm+1, n, 2)), L, t)
        assert continuant(a, c) == f*g
        assert sp.Poly(f.as_expr().subs(t, 0), L) == sp.Poly(
            sp.prod(L+a[i-1] for i in selected), L)
        assert 2*mm*mm > 8*mm-1
        print(f'FAMILY m={mm},n={n},d={2*mm}: exact membership; ambient={8*mm-1}, kappa={2*mm*mm}')
    # Check the affine spin deformation at small sizes, with symbolic parameters.
    for n in (4, 6):
        a0, c0 = spin_data(n)
        a = [u+v*ai for ai in a0]
        c = [z*ci for ci in c0]
        f = sp.prod((L+u)**2-q*q*(v*v+z*t) for q in range(1, n, 2))
        equal(continuant(a, c), f)
        assert sp.cancel((a[0]+a[-1])/2-u) == 0
        assert sp.cancel((a[0]-a[-1])/(2*(n-1))-v) == 0
        assert sp.cancel(c[0]/(n-1)-z) == 0
        print(f'AFFINE_SPIN n={n}: symbolic identity and parameter inverse passed')
    assert 31-3 == 28 < 32
    print('STRONGER_DIMENSION n=16: 3-dimensional affine spin family gives component codimension <=28<32')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
