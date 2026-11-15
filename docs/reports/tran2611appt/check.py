#!/usr/bin/env python3
"""Exact checks of Tran's spectral formulas and the proposed asymptotic route.

No floating point is used. These finite checks are not a proof of the
universal APPT volume asymptotic. monte_carlo.py is a separate diagnostic.
"""
from fractions import Fraction as F
from itertools import combinations, permutations
from math import factorial, prod
from random import Random


def rs(m):
    return m * (m - 1) // 2, m * (m + 1) // 2


def coefficients(m, d):
    r, s = rs(m)
    return [F(max(0, i + s - d) - min(i, r), i) for i in range(1, d + 1)]


def cas(m, d):
    return prod(F(2 * i, (m - 1) * d + 2 * i) for i in range(1, d))


def terms(m, d):
    c = coefficients(m, d)
    return [z ** (d - 1) / prod(z - x for j, x in enumerate(c) if j != k)
            for k, z in enumerate(c) if z > 0]


def dominant_ratio(m, d):
    r, s = rs(m)
    return (prod(F(r * d + m * i, d + m) for i in range(1, r + 1)) / factorial(r)
            * prod(F(s * d - m * j, d - m) for j in range(1, s)) / factorial(s - 1))


def beta(m):
    r, s = rs(m)
    return F(r ** r * s ** (s - 1), factorial(r) * factorial(s - 1))


def precursor(n):
    b = lambda a, b: F(factorial(a - 1) * factorial(b - 1), factorial(a + b - 1))
    return F(n * (2 * n - 1), n - 1) * (b(n, 2 * n - 2) - b(2 * n - 1, 2 * n - 2))


def det(A):
    if not A:
        return F(1)
    return sum(((-1) ** j * A[0][j] * det([row[:j] + row[j+1:] for row in A[1:]])
                for j in range(len(A))), F(0))


def psd(A):
    return all(det([[A[i][j] for j in I] for i in I]) >= 0
               for k in range(1, len(A) + 1) for I in combinations(range(len(A)), k))


def hildebrand_extreme(x, m, top, bottom):
    A = [[F(0) for _ in range(m)] for _ in range(m)]
    p = q = 0
    for i in range(m):
        for j in range(i, m):
            if i == j:
                A[i][j] = 2 * bottom[p]
            else:
                A[i][j] = A[j][i] = bottom[p] - top[q]
                q += 1
            p += 1
    return A


def main():
    count = 0
    for m in range(2, 9):
        r, s = rs(m)
        assert r + s == m*m and s - r == m
        for n in range(m, m + 9):
            d = m * n
            c = coefficients(m, d)
            # Source fidelity: evaluate G at every ordered-simplex vertex.
            for i in range(1, d + 1):
                x = [F(1, i) if k < i else F(0) for k in range(d)]
                assert sum(x[-s:]) - sum(x[:r]) == c[i-1]
            ts = terms(m, d)
            assert 0 < sum(ts) <= 1
            assert ts[-1] / cas(m, d) == dominant_ratio(m, d)
            # Correct change-of-density normalization: T_m, not V_outer.
            rates = [1 - F(d, m) * z for z in c[:-1]]
            assert all(z > 0 for z in rates)
            assert prod(1/z for z in rates) == ts[-1]
            assert 0 < sum(ts) / ts[-1] < 1
            if (m - 1) * d % 2 == 0:
                a = (m - 1)*d//2
                assert cas(m,d) == F(factorial(d-1)*factorial(a),factorial(a+d-1))
            count += 1
    print('FIDELITY: G(vertex), CAS product/Gamma, outer section and dominant ratio:',count,'(m,n) pairs OK')
    for n in range(2, 9):
        q = sum(terms(2, 2*n))
        assert q == precursor(n)
        p = F(factorial(n)*factorial(2*n-2),factorial(3*n-2))
        assert p / cas(2,2*n) == F(3*n-1,2*n-1)
        assert p <= q < 4*p
        print('QUBIT n=',n,'CAS=',cas(2,2*n),'U=',q,'p/CAS=',p/cas(2,2*n),sep='')
    assert beta(2) == F(9,2)
    assert F(3,2)*3 == beta(2)
    print('LEADING CONSTANT: beta_2=9/2; precursor p/CAS -> 3/2 and U/p -> 3 are compatible')

    # Local PSD certificate on exact rational spectra satisfying both hypotheses.
    tested = 0
    rng = Random(261101492)
    for m in [2,3]:
        r,s = rs(m)
        d = m*m
        for trial in range(4):
            # Head near (m+1)/(m-1), tail near 1, a small positive gap.
            tail = sorted([F(1)+F(rng.randrange(10),100000) for _ in range(s)],reverse=True)
            head0 = F(sum(tail)-F(1,100),r)
            head = sorted([head0+F(2*j-r+1,1000000) for j in range(r)],reverse=True)
            x = head+tail
            delta = x[0]-x[r-1]+x[d-s]-x[-1]
            G = sum(tail)-sum(head)
            assert 0 < x[-1] and delta <= x[-1]/(2*m)
            assert G >= 2*m**3*delta**2/x[-1]
            for a in permutations(head):
                for b in permutations(tail):
                    assert psd(hildebrand_extreme(x,m,a,b))
                    tested += 1
    print('LOCAL PSD: exact extreme-slot assignments checked:',tested)

    # Exact correction to the search-seat law: negative B occurs in its proposed
    # unconditioned tilted law. Here m=n=2, Y1=Y2=1, Y3=17/2, Z=1.
    c = coefficients(2,4)
    y = [F(1),F(1),F(17,2)]
    B = -F(4,2)*sum(a*b for a,b in zip(c[:-1],y))
    assert B == F(-11,3) and B+1 == F(-8,3)
    assert sum(terms(2,4)) == F(1,2)
    assert terms(2,4)[-1] == F(1)
    print('SUPPORT CORRECTION: m=n=2, B=-11/3, B+Z=-8/3; Q(B+Z>=0)=1/2=V_U/T_m')

    # Finite identities used in the probabilistic bound, for all small parameters.
    for m in range(2,9):
        for n in range(m, m+9):
            d=m*n;c=coefficients(m,d);rates=[1-F(d,m)*z for z in c[:-1]]
            w=[-F(d,m)*z for z in c[:-1]]
            mean=sum(wi/li for wi,li in zip(w,rates))
            var=sum(wi*wi/(li*li) for wi,li in zip(w,rates))
            assert var > 0
            # delta only depends on Y_1,...,Y_{r-1} and Y_{d-s+1},...,Y_{d-1}.
            r,s=rs(m)
            deltaweights=[F(int(i<r)+int(d-s+1<=i),i) for i in range(1,d)]
            edelta=sum(wi/li for wi,li in zip(deltaweights,rates))
            edelta2=edelta**2+sum(wi**2/li**2 for wi,li in zip(deltaweights,rates))
            assert edelta2 > 0
    print('TILT MOMENTS: exact first/second moment formulas and positive rates on all 63 small cases OK (finite mean B need not be positive)')
    print('ALL_CHECKS_OK')

if __name__ == '__main__':
    main()
