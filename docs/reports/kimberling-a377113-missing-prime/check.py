#!/usr/bin/env python3
"""Exact arithmetic for Kimberling's OEIS A377113 missing-prime candidate.

Source: https://oeis.org/A377113, revision 10.
All computations use Python integers. No floating point or external packages.
"""
from math import comb, isqrt

P = 167
PERIOD = 166
SOURCE_TERMS = [
    1, 3, 14, 72, 400, 2328, 13904, 84192, 513472, 3143232,
    19278464, 118359552, 727045120, 4467233280, 27452300288,
    168714381312, 1036914921472, 6372994560000, 39169586880512,
    240744913207296, 1479676193996800, 9094462273585152,
    55896907276156928,
]
EXPECTED_RESIDUES = [
    1,3,14,72,66,157,43,24,114,125,151,139,99,147,111,60,10,
    143,34,86,134,99,26,70,9,70,4,131,152,152,120,90,137,151,
    62,60,157,88,29,150,159,76,12,121,166,23,45,114,150,130,
    124,39,114,131,88,74,41,3,150,53,19,110,53,138,80,149,
    139,56,10,9,128,64,157,107,142,34,51,108,150,102,49,156,
    134,83,141,42,24,27,74,19,130,150,89,25,115,74,78,60,90,
    104,124,42,92,28,131,82,139,63,125,44,166,47,46,124,42,
    137,104,65,110,39,81,46,74,58,129,77,127,116,31,143,89,
    31,61,108,144,152,158,145,63,3,144,102,45,149,11,59,54,
    106,146,152,85,17,85,57,110,22,9,33,99,10,92,2,87,108,
    137,82,
]
REPORTED_GAP_PATTERNS = {
    2: [2],
    3: [8, 1, 4, 3, 8],
    5: [12, 12, 1, 6, 4, 1, 5, 7, 12],
    7: [3, 15, 9, 9, 15, 3, 18],
}


def transition(s):
    a, b, c, d = s
    return 3*a + 2*b + 3*c, a + 3*b + 3*d, a + 3*c + 2*d, b + c + 3*d


def multinomial_coefficients(n):
    """Expand (3+x+y)^n and reduce x^2=2, y^2=3, independently of T."""
    out = [0, 0, 0, 0]
    for j in range(n + 1):
        for k in range(n - j + 1):
            scalar = comb(n, j) * comb(n-j, k) * 3**(n-j-k)
            scalar *= 2**(j//2) * 3**(k//2)
            out[(j % 2) + 2*(k % 2)] += scalar
    return tuple(out)


def main():
    assert P > 1 and all(P % d for d in range(2, isqrt(P) + 1))
    print('prime=167; trial divisors=2..12; all remainders nonzero')
    state = (1, 0, 0, 0)
    for n, expected in enumerate(SOURCE_TERMS):
        expanded = multinomial_coefficients(n)
        assert state == expanded
        assert state[0] == expected
        state = transition(state)
    print('fidelity: all 23 published DATA terms at n=0..22 equal exact multinomial expansion and T')
    print('fidelity: a(2)=14; a(3)=72; source FORMULA/EXAMPLE indexing errors are not used')
    actual_patterns = {2: [1], 3: [2], 5: [8, 1, 4, 3, 8], 7: [12, 10, 26]}
    for prime, pattern in actual_patterns.items():
        state = (1, 0, 0, 0)
        zeros = []
        for n in range(1001):
            if state[0] == 0:
                zeros.append(n)
            state = tuple(x % prime for x in transition(state))
        gaps = [y-x for x, y in zip(zeros, zeros[1:])]
        assert len(gaps) >= 5*len(pattern)
        assert all(g == pattern[i % len(pattern)] for i, g in enumerate(gaps))
        assert gaps[:len(REPORTED_GAP_PATTERNS[prime])] != REPORTED_GAP_PATTERNS[prime]
        print(f'source erratum: p={prime}; first zero index={zeros[0]}; {len(gaps)} computed gaps through n=1000 match {pattern}, not the displayed {REPORTED_GAP_PATTERNS[prime]}')
    states = []
    state = (1, 0, 0, 0)
    for _ in range(PERIOD):
        states.append(state)
        state = tuple(x % P for x in transition(state))
    residues = [s[0] for s in states]
    assert len(EXPECTED_RESIDUES) == PERIOD
    assert residues == EXPECTED_RESIDUES
    assert state == (1, 0, 0, 0)
    assert state not in states[1:]
    assert all(residues)
    recurrence = [1, 3, 14, 72]
    for n in range(4, PERIOD + 4):
        recurrence.append((12*recurrence[-1] - 44*recurrence[-2] +
                           48*recurrence[-3] + 8*recurrence[-4]) % P)
    assert recurrence[:PERIOD] == residues
    assert recurrence[PERIOD:] == recurrence[:4]
    assert 13**2 % P == 2 and 62**2 % P == 3
    signs = [(1,1), (1,-1), (-1,1), (-1,-1)]
    bases = [(3+e*13+f*62) % P for e, f in signs]
    assert bases == [78, 121, 52, 95]
    assert all(bases)
    assert all(pow(x, PERIOD, P) == 1 for x in bases)
    inverse_four = pow(4, -1, P)
    for n, s in enumerate(states):
        for (e, f), x in zip(signs, bases):
            evaluation = (s[0] + s[1]*e*13 + s[2]*f*62 + s[3]*e*f*13*62) % P
            assert evaluation == pow(x, n, P)
    powers = [sum(pow(x, n, P) for x in bases)*inverse_four % P
              for n in range(PERIOD)]
    assert powers == residues
    print('period=166; exact full state return=(1, 0, 0, 0); no earlier full state return')
    print('zero positions in a(0)..a(165) modulo 167: []')
    print('residues=' + ','.join(map(str, residues)))
    print('independent quartic recurrence: 12,-44,48,8; all 166 residues agree')
    print('split-field roots: sqrt2=13, sqrt3=62; conjugates=78,121,52,95')
    print('split-field power sum and all four evaluations agree on all 166 states')
    print('certificate: deterministic T returns to its initial state after 166 steps; all first coordinates on that cycle are nonzero')
    print('scope: exact finite certificate for an all-index argument; no Lean proof is supplied')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
