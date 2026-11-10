#!/usr/bin/env python3
"""Exact finite checks for the spatial branch of Mayeli Remark 9.4.

Only integers and Fraction are used. A pair (u,v) represents u+v/pi^2.
No finite check is asserted to prove the asymptotic existential claim.
"""
from fractions import Fraction as Q


def floor_root(n, p):
    assert n >= 0 and p >= 1
    lo, hi = 0, 1
    while hi**p <= n:
        hi *= 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if mid**p <= n:
            lo = mid
        else:
            hi = mid
    return lo


def atan_interval(x, terms=32):
    """Alternating series enclosure for atan(x), 0<x<1."""
    s = sum(((-1)**j * x**(2*j+1) / (2*j+1) for j in range(terms)), Q(0))
    t = (-1)**terms * x**(2*terms+1) / (2*terms+1)
    return min(s, s+t), max(s, s+t)


# Machin's identity pi = 16 atan(1/5) - 4 atan(1/239).
l5, u5 = atan_interval(Q(1, 5))
l239, u239 = atan_interval(Q(1, 239))
PI_LO, PI_HI = 16*l5 - 4*u239, 16*u5 - 4*l239
assert Q(333, 106) < PI_LO < PI_HI < Q(355, 113)


def kernel(k):
    """K_{(-1/4,1/4)}(k), as constant + coefficient/pi."""
    if k == 0:
        return Q(1, 2), Q(0)
    k = abs(k)
    return Q(0), (Q((-1)**((k-1)//2), k) if k % 2 else Q(0))


def direct_defect(A):
    # Actual diagonal of T-T*T, computed from individual matrix entries.
    u, v = Q(0), Q(0)
    for m in A:
        u += Q(1, 2)
        for n in A:
            a, b = kernel(m-n)
            c, d = kernel(n-m)
            assert a*d + b*c == 0  # no 1/pi term survives here
            u -= a*c
            v -= b*d
    return u, v


def correlation_defect(A):
    # Source Proposition 4.1: #A |S| - sum_k |K(k)|^2 # (A cap (A-k)).
    A = set(A)
    span = max(A)-min(A) if A else 0
    base = min(A) if A else 0
    bits = sum(1 << (n-base) for n in A)
    u, v = Q(len(A), 2), Q(0)
    for k in range(-span, span+1):
        count = (bits & (bits >> abs(k))).bit_count()
        a, b = kernel(k)
        assert a*b == 0
        u -= count*a*a
        v -= count*b*b
    return u, v


def runs(A):
    A = set(A)
    return sum(n-1 not in A for n in A)


def check_fidelity():
    subsets = 0
    for mask in range(1 << 7):
        A = [n-3 for n in range(7) if mask >> n & 1]
        assert direct_defect(A) == correlation_defect(A)
        assert len(set(A) ^ {n-1 for n in A}) == 2*runs(A)
        subsets += 1
    for N in range(1, 32):
        # Source Proposition 9.3, one-dimensional half-band specialization.
        expected = (Q(N, 4), -2*sum((Q(N-k, k*k) for k in range(1, N, 2)), Q(0)))
        assert direct_defect(range(N)) == expected
    # The full infinite symmetric-difference identity reduces exactly to this
    # finite formula, using sum_{k odd, k in Z} 1/k^2 = pi^2/4.
    print(f"FIDELITY: {subsets} subsets, 31 interval matrices: exact trace/correlation formulas agree")


def digitization(a, R):
    """Exactly (R E_a) cap Z; E_a = union_j (1/(4(2j)^a),1/(4(2j-1)^a))."""
    A = []
    for n in range(1, int(R/4)+1):
        t = R/(4*n)
        m = floor_root(t.numerator // t.denominator, a)
        # All endpoints are excluded. Otherwise x_{m+1}<n/R<x_m.
        if Q(m**a) != t and m % 2:
            A.append(n)
    return A


def check_family():
    cases, shift_checks = 0, 0
    example = None
    for a in [1, 2, 3, 4, 7]:
        for R in [Q(31, 3), Q(64), Q(257, 2), Q(1024), Q(4096), Q(16384)]:
            A = digitization(a, R)
            B = set(A)
            t = a*R/16
            M = floor_root(t.numerator//t.denominator, a+1)
            for j in range(1, M):
                assert R*(Q(1, 4*j**a)-Q(1, 4*(j+1)**a)) >= 4
            assert runs(A) >= M//2
            # Check the upper estimate without taking any irrational powers:
            # delta <= 10 R^(1/(a+1)) |k|^(a/(a+1)).
            for k in list(range(1, 65))+[127, 257, 1024, 20000]:
                delta = len(B ^ {n-k for n in B})
                assert delta**(a+1) <= 10**(a+1)*R*k**a, (a, R, k, delta)
                shift_checks += 1
            # Evaluate trace using pair differences; do not run a quadratic
            # matrix multiplication on the larger family examples.
            u, v = correlation_defect(A)
            assert v <= 0
            pi2_defect_lower = u*PI_LO**2 + v
            assert pi2_defect_lower >= 2*runs(A), (a, R)
            if M >= 2:
                # Candidate lower constant: D*pi^2 >= (aR/16)^(1/(a+1))/4.
                assert (4*pi2_defect_lower)**(a+1) >= a*R/16
            if a == 1 and R == 64:
                example = (A, u, v, M, runs(A))
            cases += 1
    print(f"SPATIAL FAMILY: {cases} rational-scale cases, {shift_checks} shifts; gap/run/trace bounds pass")
    A, u, v, M, nr = example
    print(f"EXAMPLE: a=1, gamma=1/2, R=64, A={A}, M={M}, runs={nr}, D={u}+({v})/pi^2")
    # Tensoring by the identity on the transverse spatial box scales the
    # defect by its dimension. Check the matrix operation independently.
    A = digitization(1, Q(64))
    for N in range(1, 5):
        block_points = [(x, y) for x in A for y in range(N)]
        u = Q(len(block_points), 2)
        v = Q(0)
        for x, y in block_points:
            for z, w in block_points:
                if y == w:
                    aa, bb = kernel(x-z)
                    u -= aa*aa
                    v -= bb*bb
        p, q = direct_defect(A)
        assert (u, v) == (N*p, N*q)
    print("TENSOR: four transverse identity sizes reproduce the normalized concentration matrix")


if __name__ == '__main__':
    check_fidelity()
    check_family()
    print('ALL_CHECKS_OK')
