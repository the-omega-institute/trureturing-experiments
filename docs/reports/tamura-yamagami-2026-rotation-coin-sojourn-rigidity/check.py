#!/usr/bin/env python3
"""Exact path checks for arXiv:2609.05033v1, Sections 2--5.

No floating point is used. Symbolic identities are polynomial identities over Q;
finite balanced checks do not substitute for the all-time proof.
"""
from itertools import product
from collections import defaultdict
import sys
import sympy as S

ZERO = (0, 0, 0, 0)
ID = (1, 0, 0, 1)

def add(a, b):
    return tuple(S.expand(x + y) for x, y in zip(a, b))

def scale(c, a):
    return tuple(S.expand(c * x) for x in a)

def mul(a, b):
    return tuple(S.expand(x) for x in (
        a[0]*b[0]+a[1]*b[2], a[0]*b[1]+a[1]*b[3],
        a[2]*b[0]+a[3]*b[2], a[2]*b[1]+a[3]*b[3]))

def matrix(a):
    return S.Matrix(2, 2, a)

def split(u):
    return (u[0], u[1], 0, 0), (0, 0, u[2], u[3])

def positive(x, d):
    return int(x > 0 or (x == 0 and d == 1))

def recurrence(u, n):
    p, q = split(u)
    states = {(0, 0): ID}
    for _ in range(n):
        nxt = defaultdict(lambda: ZERO)
        for (x, k), a in states.items():
            for d, m in [(-1, p), (1, q)]:
                key = (x+d, k+positive(x, d))
                nxt[key] = add(nxt[key], mul(m, a))
        states = dict(nxt)
    return [states.get((0, k), ZERO) for k in range(n+1)]

def enumerate_words(u, n):
    # Independent implementation: inspect all complete words and their prefixes.
    p, q = split(u)
    out = [ZERO for _ in range(n+1)]
    returns = 0
    for word in product((-1, 1), repeat=n):
        if sum(word):
            continue
        prefixes = [sum(word[:j]) for j in range(n)]
        k = sum(int(x > 0 or (x == 0 and word[j] == 1))
                for j, x in enumerate(prefixes))
        a = ID
        for d in word:
            a = mul(p if d == -1 else q, a)
        out[k] = add(out[k], a)
        returns += 1
    return out, returns

def equal(a, b):
    assert matrix(a).applyfunc(S.expand) == matrix(b).applyfunc(S.expand), (a, b)

def gram(a):
    m = matrix(a)
    return (m.conjugate().T*m).applyfunc(S.expand)

def main():
    print('Python', sys.version.split()[0], 'SymPy', S.__version__)
    c, s = S.symbols('c s', real=True)
    r = (c, s, s, -c)
    g = recurrence(r, 8)
    enum, count = enumerate_words(r, 8)
    assert count == 70
    for a, b in zip(g, enum):
        equal(a, b)
    print('70 symbolic returning words: enumeration = recurrence, all nine bins')

    f = c**4-3*c*c*s*s+s**4
    gg = 3*c**4-5*c*c*s*s+s**4
    h = c**4-5*c*c*s*s+3*s**4
    j = c**4-4*c*c*s*s+2*s**4
    k = 2*c**4-4*c*c*s*s+s**4
    expected = {
        2: (-c*c*s*s*f, -c*s**3*f, c*s**3*gg, -c*c*s*s*h),
        4: (-c*c*s*s*j, -c*s**3*k, c*s**3*k, -c*c*s*s*j),
        6: (-c*c*s*s*h, -c*s**3*gg, c*s**3*f, -c*c*s*s*f)}
    for l, a in expected.items():
        equal(g[l], a)
    obstruction = (gram(g[2])+gram(g[6])-2*gram(g[4])).applyfunc(S.factor)
    scalar = 2*c**4*s**6*(c*c-s*s)**2*(c*c+s*s)
    assert (obstruction-scalar*S.eye(2)).applyfunc(S.expand) == S.zeros(2)
    print('Exact polynomial Gram certificate:', S.factor(scalar), '* I_2')

    # Fidelity: source Section 4, Hadamard n=8 matrix, weights and normalization.
    integer_h = (1, 1, 1, -1)
    gh = [scale(S.Rational(1, 16), a) for a in recurrence(integer_h, 8)]
    m = (1, 1, -1, 1)
    psi = S.Matrix([1, S.I])/S.sqrt(2)
    weights = []
    for l, a in enumerate(gh):
        equal(a, scale(S.Rational(1, 16), m) if l in (2,4,6) else ZERO)
        v = matrix(a)*psi
        weights.append(S.simplify((v.conjugate().T*v)[0]))
    assert weights == [0,0,S.Rational(1,128),0,S.Rational(1,128),0,S.Rational(1,128),0,0]
    z = sum(weights)
    assert z == S.Rational(3,128)
    assert [w/z for w in weights] == [0,0,S.Rational(1,3),0,S.Rational(1,3),0,S.Rational(1,3),0,0]
    print('Source fidelity: Gamma_8(2,4,6)=M/16; w=1/128; Z=3/128; mu=1/3')

    for mindex in range(1, 9):
        n = 4*mindex
        gh = [scale(S.Rational(1, 2**(n//2)), a) for a in recurrence(integer_h, n)]
        b = S.binomial(S.Rational(1,2), mindex)
        for l, a in enumerate(gh):
            equal(a, scale(-b/2, m) if l%2 == 0 and 2 <= l <= n-2 else ZERO)
        assert b != 0
        assert gram(scale(-b/2, m)) == b*b/2*S.eye(2)
    print('Balanced matrix coefficient and state-independent Gram weights: m=1,...,8 (n<=32)')

    rr = (S.Rational(3,5), S.Rational(4,5), S.Rational(4,5), -S.Rational(3,5))
    ss = S.diag(S.I, 1)
    u = tuple((S.diag(-1, S.I)*ss*matrix(rr)*ss.conjugate().T))
    assert matrix(u).conjugate().T*matrix(u) == S.eye(2)
    for n in range(2, 13, 2):
        gu, gr = recurrence(u,n), recurrence(rr,n)
        for a,b in zip(gu,gr):
            equal(a, tuple((-S.I)**(n//2)*ss*matrix(b)*ss.conjugate().T))
    print('Exact rational-complex phase covariance: even n=2,...,12, every bin')

    assert recurrence(ID,8) == [ZERO]*9
    flip = recurrence((0,1,1,0),8)
    for l,a in enumerate(flip):
        if l not in (0,8):
            equal(a,ZERO)
    assert gram(flip[0])+gram(flip[8]) == S.eye(2)
    print('Degeneracies: diagonal Z_8=0; anti-diagonal support {0,8}, total Gram=I_2')

    # All-time sufficiency algebra: verify the rational inverse identity modulo
    # the two exact primitive-excursion quadratics; coefficient extraction is prose.
    x,y,a,b = S.symbols('x y a b')
    bp,bm = S.Matrix([[1,-1],[0,0]]),S.Matrix([[0,0],[1,1]])
    ident = S.eye(2)
    ff = (ident-b*bp-a*bm).inv()
    minus = (ident-(b-y)*bp-(a-x)*bm).inv()
    target = ident+(y*a-x*b)/(x-y)*matrix(m)
    basis = S.groebner([a*a-(1+x)*a+x/2,b*b-(1+y)*b+y/2],a,b,x,y,domain=S.QQ)
    for entry in (ff+minus)/2-target:
        numerator = S.fraction(S.cancel(entry))[0]
        assert basis.reduce(numerator)[1] == 0
    print('Formal-series inverse identity: exact polynomial reduction by excursion quadratics')
    print('ALL_CHECKS_OK')

if __name__ == '__main__':
    main()
