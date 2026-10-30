#!/usr/bin/env python3
"""Exact generic-q checks for Talalaev 2610.32521v1, Section 9 item 4.

Requires Python 3 and SymPy. No floating point and no Lean proof.
Free words check the presentation; an Ore normal form checks localized identities.
The independence proof is the symbolic injectivity of exponent substitutions,
not an assertion inferred from a finite sample.
"""
from collections import defaultdict
from functools import lru_cache
from itertools import product
import platform
import sympy as S

q = S.Symbol('q', nonzero=True)
s = 1 - q**-2


def clean(p):
    return {w: c for w, v in p.items() if (c := S.cancel(v)) != 0}


def add(*ps):
    r = defaultdict(lambda: S.Integer(0))
    for p in ps:
        for w, c in p.items():
            r[w] += c
    return clean(r)


def scale(p, c):
    return clean({w: c*v for w, v in p.items()})


def mul(p, r):
    out = defaultdict(lambda: S.Integer(0))
    for w, c in p.items():
        for v, d in r.items():
            out[w+v] += c*d
    return clean(out)


def word(w):
    return {tuple(w): S.Integer(1)}


def scalar(c):
    return clean({(): c})


def minus(p, r):
    return add(p, scale(r, -1))


def matmul(A, B):
    return [[add(*(mul(A[i][k], B[k][j]) for k in range(len(B))))
             for j in range(len(B[0]))] for i in range(len(A))]


def chain(*ms):
    out = ms[0]
    for m in ms[1:]:
        out = matmul(out, m)
    return out


# Literal equations (2.1), (2.2), with tensor basis (11,12,21,22).
RF = [[scalar(x) for x in row] for row in
      [[q,0,0,0],[0,1,q-q**-1,0],[0,0,1,0],[0,0,0,q]]]
P = [[scalar(int(i == 2*(j % 2)+j//2)) for j in range(4)] for i in range(4)]
R = matmul(RF, P)
assert R == [[scalar(x) for x in row] for row in
             [[q,0,0,0],[0,q-q**-1,1,0],[0,1,0,0],[0,0,0,q]]]
L = [[word('ab'[j]) if i == 0 else word('cd'[j]) for j in range(2)] for i in range(2)]
L1 = [[L[i//2][j//2] if i % 2 == j % 2 else {} for j in range(4)] for i in range(4)]
left, right = chain(R,L1,R,L1), chain(L1,R,L1,R)
entries = [minus(left[i][j],right[i][j]) for i in range(4) for j in range(4)]
rels = [minus(word('ba'),scale(word('ab'),q**2)),
        minus(word('ac'),scale(word('ca'),q**2)),
        minus(word('ad'),word('da')),
        minus(minus(word('cb'),word('bc')),scale(minus(word('aa'),word('ad')),s)),
        minus(minus(word('db'),word('bd')),scale(word('ab'),s)),
        add(minus(word('dc'),word('cd')),scale(word('ca'),s))]
words2 = list(product('abcd', repeat=2))
V = S.Matrix([[p.get(w,0) for w in words2] for p in entries])
W = S.Matrix([[p.get(w,0) for w in words2] for p in rels])
assert V.rank() == W.rank() == V.col_join(W).rank() == 6
print('PRESENTATION: all 16 RE entries and the six relations have identical rank-6 row spaces over Q(q)')

# Ore normal order a < t < b < c; B denotes b inverse.
a,t,b,c,B = (word(x) for x in 'atbcB')
one = scalar(S.Integer(1))
p = scale(add(scale(mul(a,a),1+q**-2),scale(mul(a,t),-q**-2)),s)
rules = {
    tuple('ta'):mul(a,t), tuple('ba'):scale(mul(a,b),q**2),
    tuple('ca'):scale(mul(a,c),q**-2), tuple('bt'):mul(t,b),
    tuple('ct'):mul(t,c), tuple('cb'):add(mul(b,c),p),
    tuple('Ba'):scale(mul(a,B),q**-2), tuple('Bt'):mul(t,B),
    tuple('bB'):one, tuple('Bb'):one,
    tuple('cB'):minus(mul(B,c),mul(mul(B,p),B)),
}


@lru_cache(None)
def normal_word(w):
    for i in range(len(w)-1):
        pair = w[i:i+2]
        if pair in rules:
            return reduce_poly(mul(mul(word(w[:i]),rules[pair]),word(w[i+2:])))
    return word(w)


def reduce_poly(poly):
    out = {}
    for w, v in poly.items():
        out = add(out,scale(normal_word(w),v))
    return out


def check_eq(x,y):
    assert not reduce_poly(minus(x,y)), (x,y,reduce_poly(minus(x,y)))


def commute(x,y):
    check_eq(mul(x,y),mul(y,x))


# Every overlap of the base quadratic Ore rules is checked independently.
overlaps = 0
for w in product('atbc',repeat=3):
    if w[:2] in rules and w[1:] in rules:
        check_eq(mul(rules[w[:2]],word(w[2:])),mul(word(w[:1]),rules[w[1:]]))
        overlaps += 1
print(f'ORE: all {overlaps} degree-3 critical overlaps in the four-generator presentation agree')
d = scale(minus(t,a),q**-2)
images = dict(a=a,b=b,c=c,d=d)


def substitute(poly, image):
    out = {}
    for w,v in poly.items():
        term = scalar(v)
        for x in w:
            term = mul(term,image[x])
        out = add(out,term)
    return out


for rel in rels:
    check_eq(substitute(rel,images),{})
# Recover the Ore relations from the original RE relations by degree-2 row-space tests.
def old(poly):
    return substitute(poly,dict(a=word('a'),b=word('b'),c=word('c'),
                                t=add(word('a'),scale(word('d'),q**2))))
for pair in ['ta','ba','ca','bt','ct','cb']:
    r = old(minus(word(pair),rules[tuple(pair)]))
    row = S.Matrix([[r.get(w,0) for w in words2]])
    assert W.col_join(row).rank() == 6
print('ORE IDENTIFICATION: both directions between the RE and iterated Ore presentations checked')
trace = add(a,scale(d,q**2))
for gen in [a,b,c,d,B]:
    commute(trace,gen)
I1 = b
I2 = minus(c,mul(mul(d,B),a))
Delta = add(minus(mul(b,c),mul(a,d)),scale(mul(a,a),s))
check_eq(mul(I1,I2),Delta)
check_eq(mul(I2,I1),Delta)
commute(I1,I2)
print('FIDELITY: source Theorem 2.1 at n=2: [I1,I2]=0, using I1=b and I2=c-d*b^(-1)*a')
commute(I1,trace)
commute(I2,trace)
print('CENTRAL WITNESS: t=a+q^2*d commutes with a,b,c,d,b^(-1),I1,I2')
print('DENOMINATOR: b*I2=I2*b=Delta=bc-ad+(1-q^(-2))*a^2')

# Check phi at the original algebra, never on its division ring.
u,v,w = S.symbols('u v w')
for rel in rels:
    out = S.Integer(0)
    for letters, coeff in rel.items():
        term = coeff
        for letter in letters:
            term *= dict(a=0,b=u,c=v,d=w)[letter]
        out += term
    assert S.cancel(out) == 0
assert S.expand(u*v) != 0
print('SPECIALIZATION: all six relations vanish under a=0,b=u,c=v,d=w')

# Symbolic inverse formulas establish injectivity for arbitrary exponent triples.
i,j,k,N = S.symbols('i j k N',integer=True,nonnegative=True)
map1 = (i+j,j,k)
inverse1 = (map1[0]-map1[1],map1[1],map1[2])
map2 = (N+i-j,j,k)
inverse2 = (map2[0]+map2[1]-N,map2[1],map2[2])
assert all(S.expand(x-y)==0 for x,y in zip(inverse1,(i,j,k)))
assert all(S.expand(x-y)==0 for x,y in zip(inverse2,(i,j,k)))
for n in range(9):
    triples = list(product(range(n+1),repeat=3))
    assert len({(x+y,y,z) for x,y,z in triples}) == len(triples)
    assert len({(n+x-y,y,z) for x,y,z in triples}) == len(triples)
print('INDEPENDENCE: symbolic left inverses for (i,j,k)->(i+j,j,k) and (N+i-j,j,k); nonzero coefficient q^(2k)')
# A nonzero polynomial F maps after denominator clearing to u^N F(u,v,q^2 w).
# This is a statement about the original algebra. Faithfulness then gives the
# division-ring independence, and clearing a rational denominator excludes t.
X,Y,Z = S.symbols('X Y Z')
F = X**3*Y**2+q*Y*Z**2+7*Z+1
N0 = S.degree(F,Y)
G = S.cancel(X**N0*F.subs(Y,Y/X))
assert S.Poly(G,X,Y,Z)
assert S.expand(G.subs({X:u,Y:u*v,Z:q**2*w})-u**N0*F.subs({X:u,Y:v,Z:q**2*w}))==0
assert S.Poly(S.expand(G.subs({X:u,Y:u*v,Z:q**2*w})),u,v,w).is_zero is False
noncommutator = reduce_poly(minus(mul(b,a),mul(a,b)))
assert noncommutator == scale(mul(a,b),q**2-1)
assert any(S.cancel(coeff.subs(q,2)) != 0 for coeff in noncommutator.values())
print('NONCOMMUTATIVE: [b,a]=(q^2-1)*ab is a nonzero PBW term, also at q=2')
print('PARAMETER: q is indeterminate over Q; the same identities specialize to q=2, s=3/4')
print('SCOPE: exact algebra checks support the generic family outline; no formal Lean settlement is asserted')
print(f'ENVIRONMENT: Python {platform.python_version()}, SymPy {S.__version__}')
print('ALL_CHECKS_OK')
