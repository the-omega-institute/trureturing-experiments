#!/usr/bin/env python3
"""Exact algebra and finite regression checks for ANN Conjecture 11.1, item 9.

Source: arXiv:2610.34388v1, Proposition 5.17 and Conjecture 11.1.
SymPy is used only with exact integers, rationals, I, symbols and radicals.
Finite tests do not prove the all-dimensional upper bound.
"""
import itertools
import sympy as S


def zero(x):
    return S.cancel(S.expand(x)) == 0


def matrix_zero(a):
    assert all(zero(x) for x in a), a


def projection(p):
    matrix_zero(p.H - p)
    matrix_zero(p * p - p)


def FH(rho, p, q):
    f = S.expand(S.trace(rho * (p*q + q*p - 2*q*p*q)))
    h = S.expand(S.trace(rho * (p + q)))
    assert not f.has(S.Float) and not h.has(S.Float)
    assert S.im(f) == 0 and S.im(h) == 0
    return S.cancel(f), S.cancel(h)


def bound(f, h, contrast_squared):
    """Compare f/h <= s/(1+s) by exact rational comparisons and squaring.

    Here s^2=(1-sqrt(1-contrast_squared))/2 and 0<=contrast_squared<1.
    Positive-side guards make both squarings order preserving.
    """
    assert 0 <= contrast_squared < 1
    assert h > 0
    if f <= 0:
        return
    q = f/h
    assert q < 1, (f, h, contrast_squared)
    r2 = S.cancel((q/(1-q))**2)
    u = 1 - 2*r2
    assert u >= 0 and u*u >= 1-contrast_squared, (f, h, contrast_squared)


def rank_one(v):
    v = S.Matrix(v)
    return v*v.H/(v.H*v)[0]


def householder(v):
    return S.eye(len(v)) - 2*rank_one(v)


def fidelity():
    # Reproduce the exact coordinate identity used in Proposition 5.17.
    beta, wr, wi, c, sig = S.symbols('beta wr wi c sig', real=True)
    rho = S.Matrix([[beta, wr+S.I*wi], [wr-S.I*wi, 2-beta]])
    p = S.Matrix([[c*c, c*sig], [c*sig, sig*sig]])
    q = S.diag(1, 0)
    f, h = FH(rho, p, q)
    D = beta*(1+c*c) + (2-beta)*sig*sig
    assert zero(f - 2*c*sig*wr)
    assert zero(h - D - f)
    # Reproduce both asymptotics stated in that proposition.
    d = S.symbols('d', positive=True)
    sd = S.sqrt((1-S.sqrt(1-d*d))/2)
    value = sd/(1+sd)
    assert S.limit(value/d, d, 0, dir='+') == S.Rational(1, 2)
    assert S.simplify(S.limit(value, d, 1, dir='-') - (S.sqrt(2)-1)) == 0
    print('FIDELITY: Proposition 5.17 coordinate identity, slope 1/2 and limit sqrt(2)-1')


def symbolic_sharpness():
    # Put u=sqrt(t), t=1-2s^2; reducing modulo u^2+2s^2-1
    # keeps the entire calculation polynomial/rational.
    s, u, lam = S.symbols('s u lam', real=True)
    rho = S.Matrix([[u*u, 2*s*u], [2*s*u, 2-u*u]])
    p = S.Matrix([[1, u], [u, u*u]])/(1+u*u)
    q = S.diag(1, 0)
    projection(p)
    projection(q)
    f, h = FH(rho, p, q)
    assert zero(f - 4*s*u*u/(1+u*u))
    assert zero(h - 4*u*u*(1+s)/(1+u*u))
    assert zero(f*(1+s) - s*h)
    relation = u*u + 2*s*s - 1
    remainder = lambda x: S.rem(S.Poly(S.expand(x), u), S.Poly(relation, u)).as_expr()
    assert zero(remainder(rho.det() - u**4))
    cp = rho.charpoly()
    charpoly = cp.as_expr().subs(cp.gen, lam)
    assert zero(remainder(charpoly - ((lam-1)**2-4*s*s*(1-s*s))))
    # For 0<=s<1/sqrt(2): t>0, trace=2, det=t^2>0 and
    # eigenvalues 1+-delta with delta^2=4s^2(1-s^2)<1.
    print('SHARPNESS SYMBOLIC: projections, F, H and characteristic polynomial agree exactly')

    # The nontrivial exact instance from the proposed equality family.
    rho = S.Matrix([[7, 6*S.sqrt(7)], [6*S.sqrt(7), 43]])/25
    p = S.Matrix([[25, 5*S.sqrt(7)], [5*S.sqrt(7), 7]])/32
    q = S.diag(1, 0)
    projection(p)
    projection(q)
    assert sorted(rho.eigenvals()) == [S.Rational(1,25), S.Rational(49,25)]
    f, h = FH(rho, p, q)
    assert f == S.Rational(21,40) and h == S.Rational(7,5)
    assert f/h == S.Rational(3,8)
    assert p*q != q*p
    for n in range(2, 8):
        rn = S.diag(rho, S.eye(n-2))
        pn = S.diag(p, S.zeros(n-2))
        qn = S.diag(q, S.zeros(n-2))
        projection(pn)
        projection(qn)
        assert FH(rn, pn, qn) == (f,h)
        assert rn.det() == rho.det()
    # The included endpoint delta=0: positive density and quotient exactly zero.
    assert FH(S.eye(2), S.ones(2)/2, q) == (0, 2)
    print('SHARPNESS INSTANCE: delta=24/25; F=21/40; H=7/5; F/H=3/8; padding n=2..7; delta=0')


def small_qubit_grid():
    # Exhaust every ordered projection pair from this explicit finite set,
    # including genuinely complex projections, four densities per contrast.
    ps = [S.zeros(2), S.eye(2), S.diag(1,0), S.diag(0,1)]
    ps += [rank_one([1,z]) for z in [1,-1,2,S.I,-S.I,1+S.I]]
    us = [S.eye(2)] + [householder([1,z]) for z in [1,2,S.I]]
    ds = [S.Rational(0), S.Rational(3,5), S.Rational(4,5), S.Rational(24,25)]
    checked = excluded = 0
    for p in ps:
        projection(p)
    for d,u,p,q in itertools.product(ds,us,ps,ps):
        matrix_zero(u.H*u-S.eye(2))
        rho = u*S.diag(1+d,1-d)*u.H
        f,h = FH(rho,p,q)
        if p == S.zeros(2) and q == S.zeros(2):
            assert f == h == 0
            excluded += 1
        else:
            bound(f,h,d*d)
            checked += 1
    print(f'QUBIT GRID: {checked} admissible cases; {excluded} excluded zero pairs; exact bound holds')


def all_rank_masks():
    # All pairs of spectral masks for two fixed, different rational unitary
    # frames; every rank is present. The density has a third frame.
    ds = [S.Rational(0), S.Rational(3,5), S.Rational(24,25)]
    checked = excluded = 0
    for n in range(2,5):
        u = householder(list(range(1,n+1)))
        v = householder([1]+[S.I if j%2 else j+1 for j in range(1,n)])
        w = householder([n-j for j in range(n)])
        masks = list(itertools.product([0,1], repeat=n))
        ps = [u*S.diag(*m)*u.H for m in masks]
        qs = [v*S.diag(*m)*v.H for m in masks]
        for p in ps+qs:
            projection(p)
        for d in ds:
            spectrum = [1-d+2*d*S.Rational(j,n-1) for j in range(n)]
            assert min(spectrum) == 1-d and max(spectrum) == 1+d and min(spectrum)>0
            rho = w*S.diag(*spectrum)*w.H
            for p,q in itertools.product(ps,qs):
                f,h = FH(rho,p,q)
                if p == S.zeros(n) and q == S.zeros(n):
                    assert f == h == 0
                    excluded += 1
                else:
                    bound(f,h,d*d)
                    checked += 1
    print(f'ALL RANK MASKS n=2..4: {checked} admissible cases; {excluded} zero pairs; exact bound holds')


def compression_regression():
    # P,Q have known two-dimensional blocks; rho has off-diagonal blocks.
    p1 = rank_one([1,2]);q1 = S.diag(1,0)
    p2 = rank_one([1,S.I]);q2 = S.diag(0,1)
    p=S.diag(p1,p2);q=S.diag(q1,q2)
    u=householder([1,2,3,4])
    checked=0
    for d in [S.Rational(0),S.Rational(3,5),S.Rational(24,25)]:
        for middle in itertools.product([1-d,1,1+d],repeat=2):
            rho=u*S.diag(1-d,*middle,1+d)*u.H
            r1=rho[:2,:2];r2=rho[2:,2:]
            f,h=FH(rho,p,q)
            f1,h1=FH(r1,p1,q1);f2,h2=FH(r2,p2,q2)
            assert f==f1+f2 and h==h1+h2
            for rj,fj,hj in [(r1,f1,h1),(r2,f2,h2)]:
                assert rj.trace()>0 and rj.det()>0
                dj2=S.cancel(1-4*rj.det()/rj.trace()**2)
                assert 0<=dj2<=d*d
                bound(fj,hj,dj2)
            bound(f,h,d*d)
            if d>0:
                assert rho[:2,2:] != S.zeros(2)
            checked+=1
    print(f'COMPRESSION: {checked} same-density cases; block traces add; compressed contrasts do not grow')


def main():
    print('Arithmetic: SymPy '+S.__version__+'; exact integers/rationals/complex algebra/radicals only')
    fidelity()
    symbolic_sharpness()
    small_qubit_grid()
    all_rank_masks()
    compression_regression()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
