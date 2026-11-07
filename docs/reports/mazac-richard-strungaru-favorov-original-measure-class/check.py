#!/usr/bin/env python3
"""Exact finite checks for the original Favorov periodic-dipole construction.

Only integers and fractions are used. Cyclotomic arithmetic is in
Z[z]/(z^(N/2)+1), N a power of two. No floating point or numerical Fourier
transform occurs. These checks are finite evidence, not an analytic proof.
"""
from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path


def root_mul(v, exponent):
    """Multiply by z^exponent in the indicated cyclotomic quotient."""
    d = len(v)
    out = [0] * d
    for i, a in enumerate(v):
        if a:
            q, r = divmod(i + exponent, d)
            out[r] += a if q % 2 == 0 else -a
    return out


def source_gap_witness(M):
    """Reproduce Favorov Lemma 3, alpha=1/8, by an explicit exact witness.

    c(Y)=Y^(N/8+1) product_{m=-N/8}^{N/8}(Y-z^m), N=M^2.
    Time support misses [-M/8,M/8] modulo M; its finite Fourier transform
    vanishes at all m with |m|<=N/8. Scaling by an upper bound for the
    coefficient norms gives coefficients of modulus at most one.
    """
    N = M * M
    d, g = N // 2, N // 8
    one = [1] + [0] * (d - 1)
    coeff = [one]
    for m in range(-g, g + 1):
        nxt = [[0] * d for _ in range(len(coeff) + 1)]
        for j, v in enumerate(coeff):
            zm = root_mul(v, m)
            nxt[j] = [a - b for a, b in zip(nxt[j], zm)]
            nxt[j + 1] = [a + b for a, b in zip(nxt[j + 1], v)]
        coeff = nxt
    shift = g + 1
    assert shift + len(coeff) - 1 < N - g
    # Closed forbidden time gap includes both endpoints.
    assert all(g < shift + j < N - g for j, v in enumerate(coeff) if any(v))
    for m in range(-g, g + 1):
        val = [0] * d
        for j, v in enumerate(coeff):
            term = root_mul(v, -m * (j + shift))
            val = [a + b for a, b in zip(val, term)]
        assert not any(val), (M, m)
    assert coeff[-1] == one  # witness is nonzero
    bound = max(sum(abs(a) for a in v) for v in coeff)
    assert bound >= 1
    # |sum a_r z^r| <= sum |a_r| <= bound, since |z|=1.
    for v in coeff:
        assert Q(sum(abs(a) for a in v), bound) <= 1
    # delta_{kM+j/M}: increasing k by one shifts by exactly M.
    for k in range(-3, 4):
        for j in range(shift, shift + len(coeff)):
            assert Q((k + 1) * N + j, M) - Q(k * N + j, M) == M
    print(f'FIDELITY Lemma3 alpha=1/8 M={M}: time gap, Fourier gap, nonzero, periodic OK')


def lattice_envelope():
    # For j>=2, 1/(1+j^2) <= 1/(j(j-1)). The latter telescopes.
    for j in range(2, 1000):
        assert Q(1, 1 + j*j) <= Q(1, j*(j-1))
        assert Q(1, j*(j-1)) == Q(1, j-1) - Q(1, j)
    # Two half-open central intervals contribute 2, the two tails <=3.
    assert 2 + 2 * (Q(1, 2) + 1) == 5
    cases = 0
    for M in range(1, 17):
        # A closed interval of length 1 contains <=M+1 <=2M lattice points.
        assert M + 1 <= 2*M
        assert 5 * 2*M <= 12*M
        L = 24*M + 2
        for r in range(33):
            a = Q(r, 32*M)  # all test offsets in one lattice cell
            partial = sum((Q(1, 1 + (a-Q(l, M))**2) for l in range(-L, L+1)), Q(0))
            # For |l|>L and 0<=aM<=1, compare to M^2/(|l|-1)^2.
            # Integral comparison gives the exact rational tail majorant below.
            tail_bound = Q(2*M*M, L-1)
            assert partial + tail_bound < 12*M
            for z in [-10000, -7, 0, 13, 10000]:
                translated = a + Q(z, M)
                l = r - 17
                assert translated-Q(l+z, M) == a-Q(l, M)
            cases += 1
    print(f'UNIFORM_LATTICE {cases} rational offsets: certified two-sided tails and translate reindexing OK')


def dipole_checks():
    # f(y)=1/(1+y^2) is an auxiliary C^1 envelope test, not a Schwartz
    # function. (1+y^2)|f'(y)|<=1 follows from (|y|-1)^2>=0.
    f = lambda y: Q(1, 1+y*y)
    cases = 0
    for r in range(-80, 81):
        y = Q(r, 12)
        assert 2*abs(y) <= 1+y*y
        for s in range(1, 13):
            t = Q(s, 1024)
            # Exact divided-difference identity, including y=0.
            assert f(y+t)-f(y) == -t*(2*y+t)/((1+y*y)*(1+(y+t)**2))
            assert abs(f(y+t)-f(y)) <= t
            cases += 1
    # Family of bounded rational coefficients, all small finite truncations.
    for M in range(1, 9):
        for shift in [Q(-137), Q(-1,7), Q(0), Q(5,9), Q(251)]:
            for t in [Q(1,64*M), Q(1,16*M), Q(1,8*M)]:
                for K in [0,1,2,5]:
                    val = Q(0)
                    for k in range(-K,K+1):
                        for j in range(M*M):
                            c = Q((3*j+M)%7-3,3)
                            x = k*M+Q(j,M)+shift
                            val += c*(f(x+t)-f(x))
                    assert abs(val) <= 12*M*t
    print(f'DIPOLES {cases} exact divided differences; all M<=8, shifts, t, truncations OK')


def decay_checks():
    # M_n=32^n, tau_n=2^(-3(n+4)^2), n>=1. Cubic roots and
    # tau_n^(-2/3) are rational powers of two evaluated exactly.
    two_power = lambda e: Q(2**e) if e >= 0 else Q(1,2**(-e))
    direct = []
    fourier = []
    for n in range(1, 65):
        M = 32**n
        root = Q(1,2**((n+4)**2))
        tau = root**3
        amplitude = 1/(root**2)
        assert 0 < tau < Q(1,4*M)
        assert amplitude*tau == root
        direct.append(M*root)
        fourier.append(M*M*root)
        # Symbolic exponent differences, valid as polynomial identities.
        assert 5 - ((n+5)**2-(n+4)**2) == -2*n-4
        assert 10 - ((n+5)**2-(n+4)**2) == 1-2*n
    for n in range(1,64):
        assert direct[n]/direct[n-1] == two_power(-2*n-4) <= Q(1,64)
        assert fourier[n]/fourier[n-1] == two_power(1-2*n) <= Q(1,2)
    # Fourier counting step: M^{-1}Z has <=2RM+1 points in [-R,R].
    for M in range(1,17):
        for R in range(1,17):
            assert 2*R*M+1 <= 3*R*M
            assert 2*R*M*(2*R*M+1) <= 6*R*R*M*M
    print('DECAY n=1..64: original tau bounds, exact cube roots, geometric majorants OK')
    print('FOURIER all integer M,R<=16: quadratic variation coefficient bound OK')


def main():
    print('Arithmetic: integers and fractions; cyclotomic quotient arithmetic is exact')
    for M in [4,8,16,32]:
        source_gap_witness(M)
    lattice_envelope()
    dipole_checks()
    decay_checks()
    print('Scope: finite fidelity and estimate checks; no Lean proof or infinite-family certification')
    print('check.py sha256:',sha256(Path(__file__).read_bytes()).hexdigest())
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
