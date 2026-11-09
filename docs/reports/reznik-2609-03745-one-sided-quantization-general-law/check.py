#!/usr/bin/env python3
"""Exact fidelity and small-case checks; no floating-point arithmetic.

Source: Reznik, arXiv:2609.03745v1, Example 1 and equations (1),(2),(6).
The general-law assertion is NOT established by these finite checks.
"""
from fractions import Fraction as F
from itertools import combinations, product


def compositions(total, parts):
    if parts == 1:
        yield (total,)
    else:
        for x in range(total + 1):
            for rest in compositions(total - x, parts - 1):
                yield (x,) + rest


def floor_rung(x, ladder):
    return max((F(0),) + tuple(r for r in ladder if r <= x))


def uniform_gap(ladder, b=F(1)):
    ends = (F(0),) + tuple(ladder) + (b,)
    return sum(((v-u)**2 / (2*b) for u,v in zip(ends, ends[1:])), F(0))


def atomic_gap(support, weights, ladder):
    return sum((w*(x-floor_rung(x, ladder)) for x,w in zip(support, weights)), F(0))


def atomic_dp(support, weights, n):
    # Every nonzero selected rung can be moved to the first atom above it.
    points = (F(0),) + tuple(support)
    m = len(points)
    costs = {(i,j):sum((weights[k-1]*(points[k]-points[i])
                         for k in range(i+1,j)), F(0))
             for i in range(m) for j in range(i+1,m+1)}
    dp = {0:F(0)}
    candidates = [costs[0,m]]
    for _ in range(n):
        new = {}
        for j in range(1,m):
            vals = [v+costs[i,j] for i,v in dp.items() if i<j]
            if vals:
                new[j]=min(vals)
        dp=new
        candidates.extend(v+costs[i,m] for i,v in dp.items())
    return min(candidates)


def main():
    count=0
    # Source Example 1: G_n*=1/(2(n+1)), including the pinned zero cell.
    for n in range(0,6):
        target=F(1,2*(n+1))
        optimal=tuple(F(i,n+1) for i in range(1,n+1))
        assert uniform_gap(optimal)==target
        for widths in compositions(12,n+1):
            ds=tuple(F(x,12) for x in widths)
            actual=sum((x*x/2 for x in ds),F(0))
            variance=sum(((x-F(1,n+1))**2/2 for x in ds),F(0))
            assert actual-target==variance and actual>=target
            count+=1
    print(f'FIDELITY: source Example 1, n=0..5; {count} rational partitions; G_n*=1/(2(n+1)).')

    # Threshold atoms are assigned upward to their exactly affordable rung.
    assert floor_rung(F(1,2),(F(1,2),))==F(1,2)
    assert atomic_gap((F(1,2),),(F(1),),(F(1,2),))==0
    laws=0
    for m in range(1,6):
        support=tuple(F(j,m+1) for j in range(1,m+1))
        for masses in compositions(6,m):
            weights=tuple(F(x,6) for x in masses)
            for n in range(m+1):
                brute=min(atomic_gap(support,weights,ladder)
                          for k in range(n+1) for ladder in combinations(support,k))
                assert atomic_dp(support,weights,n)==brute
                if n==m:
                    assert brute==0
                laws+=1
    print(f'ATOMIC: {laws} exact law/rung-budget cases; dynamic program equals exhaustive optimum; full-support ladders have zero loss.')

    # A source-admissible mixed law: Q(B)=B/(1+B),
    # mu = 1_[0,1](B)/(1+B)^2 dB + (1/2) delta_1.
    # In quality coordinates: (1/2) U[0,1/2] + (1/2) delta_(1/2).
    # Its ac and atomic masses are each 1/2; the proposed constant is 1/8.
    b=F(1,2)
    Q=lambda B:B/(1+B)
    Qinv=lambda x:x/(1-x)
    assert Q(F(1))==b and Q(F(0))==0
    assert b+F(1,2)==1
    constant=b*b/2
    assert constant==F(1,8)
    mixed_checks=0
    for n in range(1,9):
        ladder=tuple(b*F(i,n) for i in range(1,n+1))
        bandwidth=tuple(Qinv(x) for x in ladder)
        assert tuple(Q(B) for B in bandwidth)==ladder
        assert bandwidth[-1]==1
        optimum=F(1,8*n)
        actual=uniform_gap(ladder,b)/2+(b-floor_rung(b,ladder))/2
        assert actual==optimum and n*actual==constant
        # Check the exact nonnegative square decomposition on all small rational
        # cell-width vectors (zero widths included as boundary closures).
        for widths in compositions(7,n+1):
            hs=tuple(b*F(x,7) for x in widths)
            t=hs[-1]
            gap=sum((x*x/(4*b) for x in hs),F(0))+t/2
            squares=sum(((x-(b-t)/n)**2/(4*b) for x in hs[:-1]),F(0))
            rhs=squares+t*F(n-1,2*n)+t*t*F(n+1,4*n)/b
            assert gap-optimum==rhs and rhs>=0
            mixed_checks+=1
    print(f'MIXED: {mixed_checks} exact square-decomposition cases, n=1..8; positive atom at B=1 is included at equality; G_n*=1/(8n), nG_n*=1/8.')
    print('SCOPE: these exact checks support fidelity and small cases; no universal asymptotic theorem or literature novelty is certified.')
    print('ALL_CHECKS_OK')


if __name__=='__main__':
    main()
