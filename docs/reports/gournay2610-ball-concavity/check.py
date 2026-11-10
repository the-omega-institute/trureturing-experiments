#!/usr/bin/env python3
"""Exact checks for Gournay Q6.9; standard library only, Python >= 3.12."""
from collections import defaultdict
from fractions import Fraction as F
from math import comb
import sys
assert sys.version_info >= (3, 12)


def step(word, side, mask):
    if mask == 0:
        return word
    if word and word[-1][0] == side:
        new = word[-1][1] ^ mask
        return word[:-1] + (((side, new),) if new else ())
    return word + ((side, mask),)


def bfs(k, radius):
    balls = [{()}]
    for _ in range(radius):
        old = balls[-1]
        balls.append(old | {step(w, s, 1 << i) for w in old
                            for s in range(2) for i in range(k)})
    for n, ball in enumerate(balls):
        assert all(sum(mask.bit_count() for _, mask in w) <= n for w in ball)
    return [len(b) for b in balls]


def recurrence(k, radius):
    # A nonempty normal form has a first factor (two choices), then alternates.
    a = [comb(k, i) if i <= k else 0 for i in range(radius + 1)]
    spheres = [1]
    for n in range(1, radius + 1):
        spheres.append(a[n] + sum(a[j] * spheres[n-j] for j in range(1, n+1)))
    balls, total = [], 0
    for s in spheres:
        total += s
        balls.append(total)
    # Exact truncated series identities, including the constant coefficients.
    for n in range(radius + 1):
        lhs = spheres[n] - sum(a[j] * spheres[n-j] for j in range(1, n+1))
        assert lhs == a[n]
        assert balls[n] - (balls[n-1] if n else 0) == spheres[n]
    return spheres, balls


def compositions(n):
    if n == 0:
        yield ()
    else:
        for j in range(1, n+1):
            for tail in compositions(n-j):
                yield (j,) + tail


def composition_ball(k, n):
    total = 1
    for r in range(1, n+1):
        for c in compositions(r):
            count = 2
            for j in c:
                count *= comb(k, j) if j <= k else 0
            total += count
    return total


def padd(a, b):
    return tuple((a[i] if i < len(a) else F(0)) +
                 (b[i] if i < len(b) else F(0)) for i in range(max(len(a), len(b))))


def pscale(a, c):
    return tuple(c*x for x in a)


def pmul(a, b):
    out = [F(0)] * (len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return tuple(out)


def binomial_poly(j):
    out = (F(1),)
    for t in range(j):
        out = pscale(pmul(out, (F(-t), F(1))), F(1, t+1))
    return out


def symbolic_ball(n):
    out = (F(1),)
    for r in range(1, n+1):
        for c in compositions(r):
            term = (F(2),)
            for j in c:
                term = pmul(term, binomial_poly(j))
            out = padd(out, term)
    return out


def evaluate(a, x):
    out = 0
    for c in reversed(a):
        out = out*x+c
    return out


def shifted(a, h):
    return tuple(sum(a[j]*comb(j,i)*h**(j-i) for j in range(i,len(a)))
                 for i in range(len(a)))


def distribution(k, n):
    p = {(): F(1)}
    for _ in range(n):
        q = defaultdict(F)
        for w, mass in p.items():
            q[w] += mass/2
            for side in range(2):
                for i in range(k):
                    q[step(w,side,1 << i)] += mass/F(4*k)
        p = dict(q)
    return p


def main():
    # Source fidelity: Section 5 support(P^n)=B_n and P^(2n)(e)=||P^n||_2^2.
    for n in range(5):
        p = distribution(2,n)
        assert sum(p.values()) == 1 and all(x > 0 for x in p.values())
        assert len(p) == bfs(2,n)[n]
        assert sum(x*x for x in p.values()) == distribution(2,2*n)[()]
    print('SOURCE_FIDELITY: support(P^n)=B_n and return/L2 identity, k=2, 0<=n<=4')
    for k, radius in [(1,8),(2,8),(3,6)]:
        assert bfs(k,radius) == recurrence(k,radius)[1]
        print(f'BRUTE_FORCE_NORMAL_FORMS k={k} radius={radius} balls={bfs(k,radius)}')
    for k in [2,7,8,12,13,14,20]:
        _, b = recurrence(k,14)
        for n in [4,6,8,10,12,14]:
            assert b[n] == composition_ball(k,n)
        print(f'COMPOSITION_RECURRENCE k={k} b4={b[4]} b6={b[6]} b8={b[8]}')
    expected = {
        4: ['1','7/6','11/4','-13/6','25/4'],
        6: ['1','37/30','38/15','-65/24','275/24','-541/40','1561/120'],
        8: ['1','533/420','799/336','-403/144','945/64','-10279/360','1561/32','-47293/1008','36389/1344'],
    }
    bpoly = {n:symbolic_ball(n) for n in [4,6,8]}
    for n, a in bpoly.items():
        assert a == tuple(map(F,expected[n]))
        for k in range(1,51):
            assert evaluate(a,k) == recurrence(k,n)[1][n]
        print(f'BALL_POLYNOMIAL b{n} ascending=' + ','.join(map(str,a)))
    determinant = padd(pmul(bpoly[4],bpoly[8]),pscale(pmul(bpoly[6],bpoly[6]),-1))
    assert determinant[0] == 0
    q = tuple(c*1209600 for c in determinant[1:])
    expected_q = (-37440,25104,564504,-1346864,-4517610,22928607,-23592558,
                  -11886597,-1419750,-47011,-27546,3561)
    assert q == expected_q
    qshift = shifted(q,13)
    expected_shift = (76409570764800,970162185850560,665423504997360,210417217701640,
                      39806099322672,4975759432566,429129570153,25855104969,
                      1074472938,29471504,481677,3561)
    assert qshift == expected_shift and all(c > 0 for c in qshift)
    print('IDENTITY: b4*b8-b6^2=k*Q(k)/1209600')
    print('Q_ASCENDING=' + ','.join(map(str,q)))
    print('Q_13_PLUS_T_ASCENDING=' + ','.join(map(str,qshift)))
    print('FAMILY_CERTIFICATE: all shifted coefficients positive; every integer k>=13 has positive determinant')
    for k in range(13,101):
        b = recurrence(k,8)[1]
        assert b[4]*b[8]-b[6]**2 == k*evaluate(q,k)/1209600
        assert b[4]*b[8]-b[6]**2 > 0
    b = recurrence(8,14)[1]
    assert (b[10],b[12],b[14]) == (44894921047,5480571416319,669043676869839)
    det = b[10]*b[14]-b[12]**2
    assert det == 693152381091672 and det > 0
    print(f'K8: b10={b[10]} b12={b[12]} b14={b[14]} determinant={det}')
    # Source's classical comparison, checked from normal-form syllables of C2*C3.
    a,bside = [1],[1]
    sphere = [1]
    for n in range(1,17):
        a.append(bside[n-1] if n>1 else 1)
        bside.append(2*a[n-1] if n>1 else 2)
        sphere.append(a[n]+bside[n])
    balls = [sum(sphere[:n+1]) for n in range(17)]
    for m in range(1,9):
        assert balls[2*m] == 7*2**m-6
    print('CLASSICAL_CONTROL: C2*C3 has b_(2m)=7*2^m-6 for 1<=m<=8')
    print('ALL_CHECKS_OK')

if __name__ == '__main__':
    main()
