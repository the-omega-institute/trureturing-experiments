#!/usr/bin/env python3
"""Exact checks for the OEIS A377109 defining power and prime 241.

Only Python integers are used. No external packages or network access are needed.
The output is computational evidence; the Lean refutation remains a separate task.
"""
from math import comb, isqrt

SOURCE_DATA = [
    1, 2, 9, 38, 185, 922, 4689, 23998, 123217, 633458, 3258489,
    16765718, 86273225, 443967370, 2284733313, 11757749038,
    60508271137, 311391065570, 1602499602537, 8246883961094,
    42440638964825, 218410733951098, 1123999345270833, 5784397706237854,
]
# The search seat's proposed certificate, checked rather than used as a generator.
NOMINATED_RESIDUES = [
    1,2,9,38,185,199,110,139,66,110,169,71,45,62,52,98,113,60,33,51,
    137,52,54,22,227,103,158,199,162,97,68,57,44,36,57,188,225,22,207,16,
    60,70,15,233,64,186,37,129,16,135,193,83,142,89,89,239,87,131,206,36,
    44,8,42,218,176,132,224,176,59,216,67,124,165,160,1,23,34,7,142,90,
    181,177,223,169,13,25,180,53,172,10,183,203,46,148,233,237,195,23,194,60,
    137,197,25,136,144,205,75,188,198,175,54,240,221,73,186,121,212,193,199,218,
    121,40,13,216,153,221,125,204,187,72,165,134,77,40,37,33,233,98,37,229,
    105,74,69,87,107,65,154,21,194,75,53,233,164,74,61,125,193,44,222,81,
    181,32,11,55,96,164,22,64,136,173,197,20,110,111,104,63,208,93,202,99,
    76,227,27,153,55,170,228,113,27,238,82,189,45,122,238,86,66,226,127,25,
    60,215,227,106,222,47,195,118,52,213,179,25,78,126,218,172,74,61,198,238,
    105,219,40,201,24,167,71,10,230,153,39,175,100,111,190,58,180,215,214,42,
]


def step(state):
    """Multiply by 2 + sqrt(2) + sqrt(3) in basis (1,sqrt(2),sqrt(3),sqrt(6))."""
    a, b, c, d = state
    return 2*a + 2*b + 3*c, a + 2*b + 3*d, a + 2*c + 2*d, b + c + 2*d


def multinomial_power(n):
    """Independent literal expansion: i copies sqrt(2), j copies sqrt(3)."""
    out = [0, 0, 0, 0]
    for i in range(n + 1):
        for j in range(n - i + 1):
            out[(i % 2) + 2*(j % 2)] += (
                comb(n, i) * comb(n-i, j) * 2**(n-i-j) * 2**(i//2) * 3**(j//2)
            )
    return tuple(out)


def recurrence(count, modulus=None):
    out = [1, 2, 9, 38]
    while len(out) < count:
        x = 8*out[-1] - 14*out[-2] - 8*out[-3] + 23*out[-4]
        out.append(x if modulus is None else x % modulus)
    return out[:count]


def main():
    states = []
    state = (1, 0, 0, 0)
    for n in range(len(SOURCE_DATA)):
        assert state == multinomial_power(n), (n, state)
        states.append(state)
        state = step(state)
    actual = [s[0] for s in states]
    assert actual == SOURCE_DATA
    assert recurrence(len(SOURCE_DATA)) == SOURCE_DATA
    print('FIDELITY: all 24 OEIS data terms agree with literal multinomial expansion,')
    print('          four-coordinate multiplication, and the quartic recurrence.')
    print('SOURCE_DATA_0_23:', SOURCE_DATA)

    # Explicitly expose source misprints rather than treating them as valid evidence.
    assert states[2] == (9, 4, 4, 2)
    assert states[3] == (38, 23, 21, 12)
    print('SOURCE_MISPRINT: the displayed example says exponent 3 and a(3)=9;')
    print('                 (9,4,4,2) is the exponent-2 vector; exponent 3 is', states[3])
    print('SOURCE_MISPRINT: formula initial indices a(3)=9,a(4)=38 should be a(2)=9,a(3)=38.')
    state = (1, 0, 0, 0)
    even_zero_indices = []
    for n in range(21):
        if state[0] == 0:
            even_zero_indices.append(n)
        state = tuple(x % 2 for x in step(state))
    assert even_zero_indices == list(range(1, 21, 2))
    print('SOURCE_DISCREPANCY: p=2 zero indices through 20:', even_zero_indices)
    print('                    their gap is 2, whereas the COMMENTS list says (4).')

    p, period = 241, 240
    trial_primes = [q for q in range(2, isqrt(p)+1)
                    if all(q % d != 0 for d in range(2, isqrt(q)+1))]
    assert all(p % d != 0 for d in range(2, isqrt(p)+1))
    print('PRIME:', p, 'trial primes <= isqrt(p):', trial_primes)
    r, s = 22, 56
    assert r*r % p == 2 and s*s % p == 3
    conjugates = [(2+r+s) % p, (2-r+s) % p, (2+r-s) % p, (2-r-s) % p]
    assert conjugates == [80, 36, 209, 165]
    assert sorted(conjugates) == sorted([80, 209, 36, 165])
    assert all(x != 0 and pow(x, period, p) == 1 for x in conjugates)
    print('SPLITTING: sqrt2=22, sqrt3=56; conjugates in (++,-+,+-,--) order:', conjugates)

    state = (1, 0, 0, 0)
    orbit = []
    for n in range(period):
        orbit.append(state)
        state = tuple(x % p for x in step(state))
    assert len(set(orbit)) == period
    assert state == (1, 0, 0, 0)
    residues = [v[0] for v in orbit]
    assert len(NOMINATED_RESIDUES) == period
    assert residues == NOMINATED_RESIDUES
    assert all(x != 0 for x in residues)
    assert residues == recurrence(period, p)

    # Verify every coordinate independently by the four sign evaluations.
    inv4 = pow(4, -1, p)
    inv4r = pow(4*r, -1, p)
    inv4s = pow(4*s, -1, p)
    inv4rs = pow(4*r*s, -1, p)
    for n, v in enumerate(orbit):
        pp, mp, pm, mm = (pow(x, n, p) for x in conjugates)
        recovered = ((pp+mp+pm+mm)*inv4 % p,
                     (pp-mp+pm-mm)*inv4r % p,
                     (pp+mp-pm-mm)*inv4s % p,
                     (pp-mp-pm+mm)*inv4rs % p)
        assert recovered == v, (n, v, recovered)
    print('ORBIT: length=240; distinct full states=240; return state=', state)
    print('ZERO_RESIDUES:', [n for n, x in enumerate(residues) if x == 0])
    print('RESIDUES_0_239:', residues)
    print('INDEPENDENT_CHECKS: nominated table, quartic recurrence, and all four')
    print('                    coordinates from split-conjugate powers agree.')
    print('CERTIFICATE_SCOPE: the exact zero-free full-state cycle certifies all indices')
    print('                   once recurrence-to-source identification and iteration')
    print('                   periodicity are proved; no Lean proof is supplied here.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
