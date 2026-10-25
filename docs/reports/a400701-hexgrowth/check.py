#!/usr/bin/env python3
"""Exact finite checks for the two OEIS conjectures; not an unbounded proof.

Sources: OEIS A400701 revision 30 and A374338 revision 15.
Only Python's standard library is needed. All evidence uses integers/Fraction.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import product

D = ((1, 0), (0, 1), (1, -1), (-1, 0), (0, -1), (-1, 1))
DATA5 = [5,13,23,33,47,61,77,97,117,139,165,191,219,251,283,317,355,393,433,
         477,521,567,617,667,719,775,831,889,951,1013,1077,1145,1213,1283,
         1357,1431,1507,1587,1667,1749,1835,1921,2009,2101,2193,2287,2385,2483,2583,2687,2791]
DATA2 = [4,8,14,24,34,46,62,78,96,118,140,164,192,220,250,284,318,354,394,434,
         476,522,568,616,668,720,774,832,890,950,1014,1078,1144,1214,1284,
         1356,1432,1508,1586,1668,1750,1834,1922,2010,2100,2194,2288,2384,2484,2584,2686,2792]


def q(p):
    x, y = p
    return x*x + x*y + y*y


def update(s):
    """Native simultaneous update, with retained sites and distinct neighbors."""
    votes = Counter((x+dx, y+dy) for x,y in s for dx,dy in D)
    return s | {p for p, c in votes.items() if c >= 2}


def independent_update(s):
    """Independent candidate-window scan; membership predicates, no Counter."""
    xs, ys = zip(*s)
    return s | {(x,y) for x in range(min(xs)-1, max(xs)+2)
                for y in range(min(ys)-1, max(ys)+2)
                if sum((x+dx,y+dy) in s for dx,dy in D) >= 2}


def shape(m,n):
    return {(x,y) for x in range(-n, m+n) for y in range(-n,n+1)
            if 1-n <= 2*x+y <= 2*m-3+n
            and -n <= x-y <= m-1+n and -n <= x+2*y <= m-1+n}


def closed2(n):
    correction = (0,4,2)[n % 3]
    num = 3*n*n + 5*n + correction
    assert num % 3 == 0
    return num // 3


def circle_intersections(d):
    """Full intersections of two Euclidean unit circles, in axial coordinates.

    For d=(a,b), q=a²+ab+b². The Euclidean perpendicular to d becomes
    (-(a+2b),2a+b)/sqrt(3) in axial coordinates. Intersections are d/2
    plus/minus t*(-(a+2b),2a+b), t²=(4-q)/(12q).
    For lattice distances at most 2 the only positive q are 1,3,4,
    so t is exactly 1/2, 1/6, or 0. No numerical square root is used.
    """
    a,b=d
    z=q(d)
    t={1:F(1,2), 3:F(1,6), 4:F(0)}[z]
    assert t*t == F(4-z,12*z)
    points={(F(a,2)-sgn*t*(a+2*b),F(b,2)+sgn*t*(2*a+b)) for sgn in (-1,1)}
    for p in points:
        assert q(p)==1 and q((p[0]-a,p[1]-b))==1
        assert p[0].denominator==p[1].denominator==1
    assert len(points)==(1 if z==4 else 2)
    return {(int(x),int(y)) for x,y in points}


# 4q=(2a+b)²+3b²<=16 implies |b|<=2 and |a|<=3.
NEAR = {d for d in product(range(-4,5), repeat=2) if 0 < q(d) <= 4}
assert {q(d) for d in NEAR} == {1,3,4} and len(NEAR)==18
INTERSECTIONS = {d: circle_intersections(d) for d in NEAR}
assert all(p in D and (p[0]-d[0],p[1]-d[1]) in D
           for d, ps in INTERSECTIONS.items() for p in ps)


def circle_update(s):
    out=set(s)
    for x,y in s:
        for dx,dy in NEAR:
            if (x+dx,y+dy) in s:
                out.update((x+px,y+py) for px,py in INTERSECTIONS[(dx,dy)])
    return out


def triangle_count(n):
    return sum(1 for v in range(n+1) for w in range(n+1-v) if (v-w)%3==0)


def main():
    states={m:[{(x,0) for x in range(m)}] for m in (2,5)}
    for m in (2,5):
        for n in range(1,101):
            prev=states[m][-1]
            nxt=update(prev)
            if n<=12:
                assert nxt==independent_update(prev), (m,n,'independent update')
            assert all(-n<=x<=m-1+n and -n<=y<=n for x,y in nxt)
            if n>=2:
                assert nxt==shape(m,n), (m,n,'shape')
            states[m].append(nxt)
    assert [len(s) for s in states[5][:len(DATA5)]]==DATA5
    assert [len(s) for s in states[2][1:len(DATA2)+1]]==DATA2
    print(f'FIDELITY: A400701 n=0..{len(DATA5)-1}: {len(DATA5)} published terms matched')
    print(f'FIDELITY: A374338 n=1..{len(DATA2)}: {len(DATA2)} published terms matched')
    print('NATIVE UPDATE CROSS-CHECK: 24 full-set equalities, horizons 1..12, two seeds')
    print('SHAPE: 198 full-set equalities, horizons 2..100, two seeds')
    for n in range(2,101):
        assert len(states[5][n])==len(states[2][n])+4*n+7
    print('A400701 IDENTITY: all 99 horizons n=2..100 matched')
    for n in range(1,101):
        assert len(states[2][n])==closed2(n)
    print('A374338 THREE RESIDUE FORMULAS: all 100 horizons n=1..100 matched')
    s={(0,0),(1,0)}
    for n in range(1,9):
        s=circle_update(s)
        assert s==states[2][n], (n,'circle bridge')
    print('EUCLIDEAN CIRCLE CHECK: 18 displacement cases and 8 full-set iterations matched')
    # Separate exact enumeration of the proposed square-minus-corners count.
    for n in range(2,101):
        L=2*n+1
        c=[(L-r)//3+1 for r in range(3)]
        assert sum(t*t for t in c)-2*triangle_count(n)==len(states[2][n])
    for n in range(6,101):
        assert triangle_count(n)-triangle_count(n-6)==2*n-3
    print('CONGRUENCE COUNT: 99 square-minus-corners counts; 95 triangle recurrences matched')
    print(f'AT n=100: A400701={len(states[5][100])}; A374338={len(states[2][100])}; difference=407')
    print('ALL_CHECKS_OK')


if __name__=='__main__':
    main()
