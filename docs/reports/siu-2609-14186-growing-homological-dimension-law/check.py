#!/usr/bin/env python3
"""Exact fidelity and finite algebra checks; no simulation or floating-point evidence."""
from fractions import Fraction as F
from itertools import combinations
from math import comb
import sympy as s


def subsets(D, m):
    return list(map(frozenset, combinations(range(1, D + 1), m)))


def even_bars(D, m):
    return [(E, i, 2 ** (m - i)) for E in subsets(D, m)
            for i in sorted(E) if set(range(1, i + 1)) <= E]


def average(xs):
    return sum(xs, F(0)) / len(xs)


def mean_variance(xs):
    mu = average(xs)
    return mu, average([(x - mu) ** 2 for x in xs])


def main():
    fidelity = counts = moments = mixtures = tails = 0
    # Source Theorem 3.1: q=0 has weights 2^(D-i), total 2^D-1;
    # its q=D-2 specialization has 2D-1 bars, as does Corollary 3.2 q=0.
    for D in range(2, 11):
        bs = even_bars(D, D)
        assert [(i, a) for E, i, a in bs] == [(i, 2 ** (D-i)) for i in range(1, D+1)]
        assert sum(a for E, i, a in bs) == 2**D-1
        assert sum(a for E, i, a in even_bars(D, 2)) == 2*D-1
        fidelity += 3
        for q in range(D-1):
            neg = comb(D-1, q+1)
            pos = comb(D-1, D-q-1)
            assert neg + pos == comb(D, q+1)
            assert F(neg, neg+pos) == F(D-q-1, D)
            assert F(pos, neg+pos) == F(q+1, D)
            # Theorem 3.3, independently enumerated subsets omitting 1.
            assert len(list(combinations(range(2,D+1), q+1))) == neg
            assert len(list(combinations(range(2,D+1), D-q-1))) == pos
            mixtures += 1
        for m in range(2, D+1):
            bs = even_bars(D, m)
            A = {i:sum(a for E,j,a in bs if j == i) for i in range(1,m+1)}
            for i in range(1,m+1):
                assert A[i] == 2**(m-i) * comb(D-i,m-i)
                ratio = F(1,2**(i-1))
                for j in range(1,i):
                    ratio *= F(m-j,D-j)
                assert F(A[i],A[1]) == ratio <= F(1,2**(i-1))
                counts += 1
            for L in range(1,m+1):
                assert F(sum(a for i,a in A.items() if i > L),sum(A.values())) <= F(1,2**(L-1))
                tails += 1
            # Several bounded positive triangular-array rows, with exact rational entries.
            rows = [[F(j,D) + 1 for j in range(1,D+1)],
                    [F(1+(j % 3),2) for j in range(1,D+1)],
                    [F(3) for j in range(D)]]
            for w in rows:
                I_D = average(w)
                for i in range(1,m+1):
                    J = list(combinations(range(i,D),m-i))
                    z = [(sum(w[:i]) + sum((w[j] for j in T),F(0)))/m for T in J]
                    mu,var = mean_variance(z)
                    bias = F(D-m,m*(D-i))*(sum(w[:i])-i*I_D) if D>i else F(0)
                    assert mu-I_D == bias
                    N,h = D-i,m-i
                    vpop = mean_variance(w[i:])[1] if N else F(0)
                    expected = F(h*(N-h),N-1)*vpop/(m*m) if N>1 else F(0)
                    assert var == expected
                    assert abs(bias) <= F(i,m)*(max(w)-min(w))
                    assert var <= max(w)**2/m
                    moments += 1
    # Source equation (3.1) and the singleton observation, including odd real powers.
    # Choose lambda_j=u_j^(k-2), so source inverse powers are exactly u_j^(-2).
    for k in range(3,9):
        for u in range(1,6):
            lam = s.Integer(u)**(k-2)
            value = s.Rational(1,k)*(lam**(-s.Rational(2,k-2)))**(-s.Rational(k-2,2))
            assert s.simplify(value-lam/k) == 0
            fidelity += 1
    # Exact symbolic leading lifetime coefficient for every parity represented by k=3..10.
    # This checks the derivative and rescaling algebra, not the infinite-dimensional theorem.
    t,w = s.symbols('t w',positive=True)
    for k in range(3,11):
        r=s.Rational(k-2,2)
        nu=t**(-r)/k
        assert s.simplify(-s.diff(nu,t)*w-r*w*t**(-r-1)/k)==0
        h=s.symbols('h',positive=True)
        assert s.simplify(s.limit(((1-h*w)**(-r)-1)/(k*h),h,0)-r*w/k)==0
        # theta=1/2 gives precisely the nominated factor 2^(k/2).
        assert s.simplify(s.Rational(1,2)**(-r-1)-2**s.Rational(k,2))==0
    # Constant-array k=4 exact proportional-regime sanity check at D=200,400.
    # Finite offsets q+2 and q+1 are retained: no premature substitution of theta*D.
    for D in (200,400):
        q=D//2
        for m in (D-q,q+2,q+1,D-q-1):
            lower=F(D,4*m)
            life_removed=F(D*D,4*m*(m-1))
            life_added=F(D*D,4*m*(m+1))
            assert lower>0 and life_removed>0 and life_added>0
            assert abs(lower-F(1,2)) < F(1,20)
            assert abs(life_removed-1) < F(1,10)
            assert abs(life_added-1) < F(1,10)
    print('SOURCE_FIDELITY_IDENTITIES',fidelity)
    print('EVEN_INDEX_COUNTS_AND_RATIOS',counts)
    print('EXACT_CONDITIONAL_MEAN_VARIANCE',moments)
    print('ODD_MIXTURE_CASES',mixtures)
    print('GEOMETRIC_TAIL_CASES',tails)
    print('LIFETIME_DERIVATIVE_AND_HALF_DIMENSION k=3..10 OK')
    print('CONSTANT_ARRAY_PROPORTIONAL_CHECK D=200,400 k=4 OK')
    print('Scope: finite identities and derivative coefficients only; no proof of Wishart limits.')
    print('ALL_CHECKS_OK')


if __name__=='__main__':
    main()
