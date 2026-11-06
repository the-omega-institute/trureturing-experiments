#!/usr/bin/env python3
"""Exact finite checks for the semilinear replacement in Keller--Riedl Q9.3.

These checks are finite evidence. The unbounded wreath-family argument is
mathematical, and is not certified by this program. Standard library only.
"""
from functools import lru_cache
from itertools import permutations, product
from random import Random


def mul(x, y):
    a, b = x % 3, x // 3
    c, d = y % 3, y // 3
    return (a*c-b*d) % 3 + 3*((a*d+b*c) % 3)


def power(x, n):
    r = 1
    for _ in range(n):
        r = mul(r, x)
    return r


Z = tuple(power(4, i) for i in range(8))  # zeta = 1+t, t^2=-1
LOG = {z: i for i, z in enumerate(Z)}
assert len(LOG) == 8 and set(LOG) == set(range(1, 9))
assert power(4, 8) == 1 and power(4, 4) != 1
assert all(power(x, 9) == x for x in range(9))
print('FIELD: F_3[t]/(t^2+1), zeta order=8, Frobenius involutive')


# Transform (d,p,f): v at i moves to p[i] as zeta^d[i]*Frob^f[i](v[i]).
def compose(a, b):  # first a, then b (the paper uses right actions)
    d, p, f = a
    e, q, h = b
    return (tuple((e[p[i]] + (3 if h[p[i]] else 1)*d[i]) % 8 for i in range(8)),
            tuple(q[p[i]] for i in range(8)),
            tuple((f[i]+h[p[i]]) % 2 for i in range(8)))


def act(a, v):
    d, p, f = a
    out = [0]*8
    for i in range(8):
        out[p[i]] = mul(Z[d[i]], power(v[i], 3 if f[i] else 1))
    return tuple(out)


I = ((0,)*8, tuple(range(8)), (0,)*8)
THETA = ((0,)*8, tuple(3*i % 8 for i in range(8)), (1,)*8)
SINGER = ((1,)+(0,)*7, tuple((i+1) % 8 for i in range(8)), (0,)*8)
MGEN = []
for i in range(7):
    d = [0]*8
    d[i], d[7] = 1, 7
    MGEN.append((tuple(d), I[1], I[2]))


def in_sigma(a):
    d, p, f = a
    r = p[0]
    return f == (0,)*8 and p == tuple((i+r) % 8 for i in range(8)) and sum(d) % 8 == r


def in_replacement(a):
    d, p, f = a
    r, b = p[0], f[0]
    return (f == (b,)*8 and p == tuple(((3 if b else 1)*i+r) % 8 for i in range(8))
            and sum(d) % 8 == r and r % 2 == b)


assert compose(THETA, THETA) == I
for g in MGEN+[SINGER]:
    assert in_sigma(compose(compose(THETA, g), THETA))
assert not in_sigma(THETA)
SGEN = MGEN+[compose(SINGER, SINGER), compose(SINGER, THETA)]
assert all(in_replacement(g) for g in SGEN)
assert in_sigma(SINGER)
s = I
for _ in range(8):
    s = compose(s, SINGER)
assert s == ((1,)*8, I[1], I[2])
print('FIDELITY Lemma 5.3: Sigma normal form sum(d)=shift, |Sigma|=8^8=16777216')
print('FIDELITY Section 5.2: (zeta_0*g_3)^8=zeta*I, Singer element order=64')
print('NORMALIZER: all 8 Sigma generator conjugates are in Sigma; theta^2=1')


rng = Random(9302610)
def random_s():
    r = rng.randrange(8)
    d = [rng.randrange(8) for _ in range(7)]
    d.append((r-sum(d)) % 8)
    b = r % 2
    return (tuple(d), tuple(((3 if b else 1)*i+r) % 8 for i in range(8)), (b,)*8)


basis = [tuple(a if i == j else 0 for i in range(8)) for j in range(8) for a in (1, 3)]
composition_checks = 0
for _ in range(2048):
    a, b = random_s(), random_s()
    c = compose(a, b)
    assert in_replacement(c)
    for v in basis:
        assert act(c, v) == act(b, act(a, v))
        composition_checks += 1
print('REPLACEMENT: 2048 exact closure checks;', composition_checks, 'F_3-basis composition checks')


# Recover the unique transporter from each of the two representative types.
def transporter(target, twisted):
    zeros = [i for i, v in enumerate(target) if v == 0]
    assert len(zeros) <= 1
    r = zeros[0] if zeros else sum(LOG[v] for v in target) % 8
    b = r % 2 if twisted else 0
    p = tuple(((3 if b else 1)*i+r) % 8 for i in range(8))
    d = [LOG[target[p[i]]] if target[p[i]] else 0 for i in range(8)]
    if zeros:
        d[0] = (r-sum(d)) % 8
    ans = (tuple(d), p, (b,)*8)
    assert (in_replacement(ans) if twisted else in_sigma(ans))
    source = ((0,)+(1,)*7) if zeros else (1,)*8
    assert act(ans, source) == tuple(target)
    return ans


for _ in range(2048):
    target = [rng.randrange(1, 9) for _ in range(8)]
    for zero in (None, rng.randrange(8)):
        w = target.copy()
        if zero is not None:
            w[zero] = 0
        transporter(w, False)
        transporter(w, True)
# Stabilizers at representatives: enumerate all 8 permutation choices and the
# only possible diagonal at the unconstrained zero coordinate.
for twisted in (False, True):
    counts = []
    for zero in (False, True):
        v = (0,)+(1,)*7 if zero else (1,)*8
        count = 0
        for r in range(8):
            b = r % 2 if twisted else 0
            p = tuple(((3 if b else 1)*i+r) % 8 for i in range(8))
            for d0 in range(8):
                a = ((d0,)+(0,)*7, p, (b,)*8)
                if (in_replacement(a) if twisted else in_sigma(a)) and act(a,v) == v:
                    count += 1
        counts.append(count)
    assert counts == [1,1]
print('LOCAL ORBITS: Sigma and S representative stabilizers=[1,1]; 8192 transporter checks')
assert 8**8 == 8*(8**7)


# A diagonal coefficient cancels in commutation with each diagonal generator.
# Enumerate every local permutation/global-Frobenius pair, not just W_3.
centralizer_pairs = []
for p in permutations(range(8)):
    for b in (0,1):
        if all(all((3 if b else 1)*g[0][i] % 8 == g[0][p[i]] for i in range(8))
               for g in MGEN):
            centralizer_pairs.append((p,b))
assert centralizer_pairs == [(I[1],0)]
scalars = []
for a in range(8):
    scalar = ((a,)*8,I[1],I[2])
    if all(compose(scalar,g)==compose(g,scalar) for g in SGEN):
        assert in_replacement(scalar)
        scalars.append(a)
assert scalars == [0,4]
projection = {0}
while True:
    more = projection | {g[1][i] for g in SGEN for i in projection}
    if more == projection:
        break
    projection = more
assert len(projection) == 8
# A semilinear permutation twists diagonal coefficients. Transitivity does
# NOT force a central diagonal to be a constant scalar. The K-linear shift
# by two forces separate constants on even and odd labels; check all 64.
center_diagonals = []
for a,b in product(range(8),repeat=2):
    d = tuple(a if i % 2 == 0 else b for i in range(8))
    z = (d,I[1],I[2])
    if all(compose(z,g) == compose(g,z) for g in SGEN):
        assert in_replacement(z)
        center_diagonals.append((a,b))
assert center_diagonals == [(a,3*a % 8) for a in range(8)]
z = (tuple(1 if i % 2 == 0 else 3 for i in range(8)),I[1],I[2])
zpower = I
for e in range(1,9):
    zpower = compose(zpower,z)
    assert (zpower == I) == (e == 8)
assert all(act(compose(z,g),v) == act(compose(g,z),v) for g in SGEN for v in basis)
print('CENTER: all 80640 permutation/Frobenius pairs and all 64 remaining diagonal candidates checked')
print('CENTER RESULT: Z(S) has 8 elements, exponent pairs (even,odd)='+str(center_diagonals))
print('CENTER RESULT: diag(zeta on even labels, zeta^3 on odd labels) is central of order 8')
print('EXPECTED-ROUTE CHECK FAILED: the proposed |Z(S)|=2 and center-order nonconjugacy obstruction are false')


# Reproduce Theorem B/4.5 by enumerating ALL binary-tree subset orbit shapes,
# rather than sampling vectors. Equal child shapes contribute one symmetry bit.
@lru_cache(None)
def shapes(k):
    if k == 0:
        return (0,1)
    prev = shapes(k-1)
    return tuple((a,b) for i,a in enumerate(prev) for b in prev[i:])


@lru_cache(None)
def shape_data(s):
    if isinstance(s,int):
        return (s,0)
    a,b = map(shape_data,s)
    return (a[0]+b[0],a[1]+b[1]+int(s[0]==s[1]))


for n in range(1,6):
    for k in range(6):
        costs = [n*size+sym for size,sym in map(shape_data,shapes(k))]
        a = min(costs)
        exponent = n*2**k+2**k-1-a
        expected = n*2**k if k<=n else (n*2**(n+1)+1 if k==n+1 else n*2**k+3*2**(k-n-2)-1)
        number = 1 if k<n or k==n+1 else (2 if k==n else 3)
        assert exponent == expected and costs.count(a) == number
print('FIDELITY Theorem B/4.5: all tree subset shapes, n=1..5, k=0..5, exponents and orbit counts agree')


def wreath(g,h,b):
    out = [0]*(2*len(g))
    for i in range(len(g)):
        out[2*i] = 2*g[i]+b
        out[2*i+1] = 2*h[i]+(1-b)
    return tuple(out)


def ysets(j):
    if j == 2:
        return ({0},{0,1},{0,1,3})
    a,b,c = ysets(j-1)
    join = lambda x,y: {2*i for i in x}|{2*i+1 for i in y}
    return (join(a,b),join(a,c),join(b,c))


def stabilizers(group,j):
    return [sum({g[i] for i in y}==y for g in group) for y in ysets(j)]


E2 = [tuple(i^c for i in range(4)) for c in range(4)]
R3 = [wreath(g,h,b) for g,h,b in product(E2,E2,(0,1))]
W2 = [p for p in permutations(range(4)) if all((i%2==j%2)==(p[i]%2==p[j]%2) for i,j in product(range(4),repeat=2))]
def chars(g):
    p0 = g[0]%2
    p1 = sum(g[i]>g[j] for i in range(4) for j in range(i+1,4))%2
    return (p0^p1,p0)
V3 = [wreath(g,h,b) for g,h,b in product(W2,W2,(0,1)) if chars(g)==chars(h)]
R4 = [wreath(g,h,b) for g,h,b in product(V3,V3,(0,1))]
for j,group in [(3,R3),(4,R4)]:
    assert len(set(group)) == 2**(3*2**(j-2)-1)
    assert stabilizers(group,j) == [2,1,2]
    assert {g[0] for g in group} == set(range(2**j))
    print('FIDELITY Corollary 5.11: j='+str(j)+', |R|='+str(len(group))+', stabilizers=[2,1,2], transitive')
# Integer identities needed in the symbolic family; this finite range is not
# a proof of their universal validity.
for j in range(3,65):
    m,k = j+5,j+3
    e = 24*2**j+3*2**(j-2)-1
    assert e == 99*2**(m-7)-1 == 3*2**k+3*2**(k-5)-1
assert 2**(8-1)==128
print('FAMILY INGREDIENTS: order exponents checked for j=3..64; first instance dim_F3=128, order=2^197')
print('LIMIT: no settlement of Question 9.3; the proposed nonconjugacy proof fails its exact center check')
print('ALL_CHECKS_OK')
