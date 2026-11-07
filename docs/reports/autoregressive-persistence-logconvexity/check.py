#!/usr/bin/env python3
"""Exact finite checks for arXiv:2604.05670v1, Remark 3(a).

No floating point, sampling, numerical quadrature, or Lean proof is used.
Fidelity: integrate killed Laplace kernels by polynomial coefficients.
Candidate: classify every band box by exact interval propagation, independently
of the forbidden-pair recognizer and of the weighted recurrence.
"""
from fractions import Fraction as F
from itertools import product
from math import comb, factorial


def laplace_fidelity(limit=12):
    # For theta=1 and density exp(-abs(x))/2, the killed n-step density on
    # z>=0 is exp(-z)*P_n(z). Split the convolution integral at x=z.
    poly = [F(1, 2)]
    surv = [F(1)]
    for n in range(1, limit + 2):
        mass = sum((c * factorial(k) for k, c in enumerate(poly)), F(0))
        surv.append(mass)
        assert mass == F(comb(2*n, n), 4**n), (n, mass)
        # P_{n+1}(z)=1/2*[integral_0^z P_n(x) dx +
        # exp(2z)*integral_z^infty exp(-2x)*P_n(x) dx].
        nxt = [F(0)] * (len(poly) + 1)
        for k, c in enumerate(poly):
            nxt[k + 1] += c / (2 * (k + 1))
            for j in range(k + 1):
                nxt[j] += c * F(factorial(k), factorial(j) * 2**(k-j+2))
        poly = nxt
    for n in range(limit + 1):
        catalan = comb(2*n, n) // (n + 1)
        assert surv[n] - surv[n+1] == F(catalan, 2**(2*n+1))
        assert surv[n] * surv[n+2] >= surv[n+1]**2 if n < limit else True
    print('FIDELITY: exact killed-Laplace integrals reproduce p_n=binom(2n,n)/4^n')
    print('FIDELITY: q_n=C_n/2^(2n+1), source Section 4.2, through n=12')


def language(word):
    return not word or (word[0] and all(word[j-1] or word[j] for j in range(1, len(word))))


def box_survives(theta, positive, negative, word):
    lo = hi = F(0)
    for sign in word:
        a, b = positive if sign else negative
        lo, hi = theta * lo + a, theta * hi + b
        if hi < 0:
            return False
        assert lo >= 0, ('ambiguous band box', theta, word, lo, hi)
    return True


def band_word_counts(theta, positive, negative, limit=12):
    a, A = positive
    negB, negb = negative
    b, B = -negb, -negB
    assert 0 < theta < 1 and 0 < b <= B and 0 < a <= A
    M = A / (1 - theta)
    margin_one = theta*a-B
    margin_two = theta**2*M-(1+theta)*b
    assert margin_one > 0 and margin_two < 0
    counts = []
    for n in range(limit+1):
        row = [0]*(n+1)
        for word in product((False, True), repeat=n):
            actual = box_survives(theta, positive, negative, word)
            assert actual == language(word), (theta, word)
            if actual:
                row[sum(word)] += 1
        counts.append(row)
    return counts, margin_one, margin_two


def weighted_counts(counts, u):
    return [sum((F(c)*u**k*(1-u)**(n-k) for k, c in enumerate(row)), F(0))
            for n, row in enumerate(counts)]


def recurrence(u, limit=66):
    p, v = [F(1), u], u*(1-u)
    for _ in range(2, limit+1):
        p.append(u*p[-1]+v*p[-2])
    return p


def mmul(A, B):
    return [[sum((A[i][k]*B[k][j] for k in range(2)), F(0)) for j in range(2)] for i in range(2)]


def candidate_checks():
    us = sorted({F(a,b) for b in range(2,10) for a in range(1,b)})
    for m in (2,3,4,5):
        theta, c, width = F(1,m), F(m*m+m-1,2), F(1,10)
        positive, negative = (c-width,c+width), (-1-width,-1+width)
        counts, lower, upper = band_word_counts(theta, positive, negative)
        for u in us:
            p = recurrence(u)
            assert weighted_counts(counts,u) == p[:len(counts)]
            v = u*(1-u)
            for n in range(65):
                assert p[n]*p[n+2]-p[n+1]**2 == (-1)**n*v**(n+1)
            for n in range(33):
                for k in range(33):
                    assert p[n+k] >= p[n]*p[k]  # Source equation (12).
            # Exact independent killed-matrix calculation; start in state '-'.
            K = [[u,1-u],[u,F(0)]]
            P = [[F(1),F(0)],[F(0),F(1)]]
            for n in range(17):
                assert sum(P[1]) == p[n]
                P = mmul(P,K)
        print(f'BANDS: theta={theta}, positive={positive}, negative={negative}')
        print(f'  one-negative lower margin={lower}; two-negative upper margin={upper}')
    print(f'FAMILY: all sign words of lengths 0..12 on four band pairs; {len(us)} rational weights each')
    print('FAMILY: adjacent-minor identity through n=64; killed-matrix row sums through n=16')
    print('SOURCE CROSS-CHECK: p_(n+m)>=p_n*p_m for 0<=n,m<=32 on these families')
    u, v = F(2,7), F(10,49)
    p = recurrence(u)
    assert p[1:4] == [F(2,7),F(2,7),F(48,343)]
    D1 = p[1]*p[3]-p[2]**2
    assert D1 == -F(100,2401)
    assert u*F(5,2)+(1-u)*(-1) == 0
    print(f'WITNESS: p1={p[1]}, p2={p[2]}, p3={p[3]}, D1={D1}; component-mean mixture=0')
    eps = F(1,1000)
    assert eps < F(25,4802)
    bound = D1+8*eps
    assert bound == -F(10099,300125) < 0
    assert 8*F(25,4802) == -D1
    print(f'ROBUSTNESS ARITHMETIC: D1+8/1000={bound}<0 (coupling bound is a paper argument)')
    # Infinite rational parameter subfamily: c=(m^2+m-1)/2, width=1/10.
    # These two expressions are the margins after substitution theta=1/m.
    for m in range(2,101):
        t, c, w = F(1,m), F(m*m+m-1,2), F(1,10)
        assert t*(c-w)-(1+w) == F(5*m*m-6*m-6,10*m) > 0
        assert t*t*(c+w)/(1-t)-(1+t)*(1-w) == -F(4*m*m-5*m-5,10*m*(m-1)) < 0
    print('PARAMETERS: rational margin formulas checked for 2<=m<=100; no unbounded inference from enumeration')


if __name__ == '__main__':
    laplace_fidelity()
    candidate_checks()
    print('ALL_CHECKS_OK')
