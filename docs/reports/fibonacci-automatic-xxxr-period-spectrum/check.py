#!/usr/bin/env python3
"""Exact, bounded checks for the first assertion of DMRSS Conjecture 11.

No floating point and no external packages. This is evidence for finite cases,
not a proof of the unbounded assertion. See README.md for the finite cutoff.
"""
from itertools import product


def complement(w):
    return w.translate(str.maketrans('01', '10'))


def fib(n):
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a


def primitive_root(w):
    for d in range(1, len(w) + 1):
        if len(w) % d == 0 and w == w[:d] * (len(w) // d):
            return w[:d]
    raise AssertionError('nonempty word required')


def suffix_ok(w):
    for k in range(1, len(w) // 3 + 1):
        a, b, c = w[-3*k:-2*k], w[-2*k:-k], w[-k:]
        if a == b and c == a[::-1]:
            return False
    return True


def periodic_witness(w):
    """Return (start, root length, root) witnessing xx x^R, or None.

    It suffices to examine starts < least period and root lengths <= it.
    """
    w = primitive_root(w)
    n = len(w)
    t = w * 4
    for k in range(1, n + 1):
        for i in range(n):
            a = t[i:i+k]
            if a == t[i+k:i+2*k] and a[::-1] == t[i+2*k:i+3*k]:
                return i, k, a
    return None


def direct_witness(w):
    """Independent modular-index implementation, used on every word <= 9."""
    n = len(w)
    for k in range(1, n + 1):
        for i in range(n):
            if all(w[(i+j) % n] == w[(i+k+j) % n]
                   and w[(i+j) % n] == w[(i+3*k-1-j) % n]
                   for j in range(k)):
                return True
    return False


def h(w):
    return ''.join(complement({'S': '00100', 'L': '00100100'}[a])
                   if i % 2 else {'S': '00100', 'L': '00100100'}[a]
                   for i, a in enumerate(w))


def phi(w):
    return ''.join({'S': 'SL', 'L': 'SLL'}[a] for a in w)


def theta(w):
    t = phi(phi(phi(w)))
    assert t.startswith('SLSL')
    return t[4:] + 'SLSL'


def rote_prefix(n):
    w = '0'
    while len(w) < n:
        w = ''.join('01' if a == '0' else '0' for a in w)
    return ''.join(('0' if i % 2 == 0 else '1') * (2 if a == '0' else 1)
                   for i, a in enumerate(w))[:n]


def suitable(a, b):
    ha, hb = h(a), h(b)
    assert len(a) % 2 == len(b) % 2 == 1
    assert len(ha) % 2 == 1 and len(hb) % 2 == 0
    ell = hb[:len(hb)//2]
    assert hb == ell + ell[::-1]
    assert ha.startswith(ell)
    mu = ha[len(ell):]
    assert mu and ell and ha == mu[::-1] + ell[::-1]
    pref = ell + mu + complement(mu[::-1])
    assert hb.startswith(pref) and len(hb) > len(pref)


def cyclic_contains(u, v):
    return v in u * ((len(u) + len(v) - 2) // len(u) + 1)


def forbidden_patterns():
    b = []
    for a in 'SL':
        for t in ['L', 'SS', 'SL']:
            b.append(a+'SSSL'+t)
            b.append('LSSL'+t)
        for c in 'SL':
            b += [a+'LLLLL'+c, a+'LSLLL'+c,
                  phi(a+'SS'+c), phi(phi(a+'LSS'+c)),
                  phi(phi(a+'SSSSS'+c))]
        for t in ['L', 'SS', 'SL']:
            b.append(phi(a+'LLLSL'+t))
        b.append(phi(phi('LLL'+a)))
    return set(b)


def cyclic_phi_preimage(u):
    i = u.index('S')
    v = u[i:] + u[:i]
    pieces = ['S'+a for a in v[1:].split('S')]
    assert all(a in ['SL', 'SLL'] for a in pieces), (u, pieces)
    return ''.join('S' if a == 'SL' else 'L' for a in pieces)


def check_cyclic_descent():
    b = forbidden_patterns()
    survivors = 0
    for n in range(1, 17):
        for letters in product('SL', repeat=n):
            u = ''.join(letters)
            if any(cyclic_contains(u, v) for v in b):
                continue
            survivors += 1
            if set(u) == {'S'}:
                continue
            v = cyclic_phi_preimage(u)
            if set(v) == {'L'}:
                continue
            v = cyclic_phi_preimage(v)
            v = cyclic_phi_preimage(v)
            assert cyclic_contains(u, theta(v)), (u, v)
            assert len(u) >= 13*len(v)
    print(f'Cyclic descent helper: all 131070 auxiliary words of length <= 16; {survivors} avoid the forbidden factors and satisfy a terminal case or three parsings')


def main():
    # Independent implementations must agree before larger checks are trusted.
    for n in range(1, 10):
        for a in product('01', repeat=n):
            w = ''.join(a)
            assert (periodic_witness(w) is not None) == direct_witness(w), w
    print('Independent modular/string periodic check: all 1022 nonempty words of length <= 9 agree')

    # Currie--Rampersad, arXiv:1502.07014v1, Table 1 (all entries n <= 48).
    source_counts = [2,4,6,10,16,24,34,48,62,80,100,124,148,178,210,244,
                     282,324,372,426,488,556,630,712,804,908,1024,1152,
                     1296,1454,1626,1814,2018,2244,2490,2756,3044,3354,
                     3690,4050,4438,4856,5300,5772,6272,6800,7370,7966]
    words = ['']
    spectrum = []
    for n, expected in enumerate(source_counts, 1):
        words = [w+a for w in words for a in '01' if suffix_ok(w+a)]
        assert len(words) == expected, (n, len(words), expected)
        cycles = [w for w in words if primitive_root(w) == w
                  and periodic_witness(w) is None]
        if cycles:
            spectrum.append(n)
            assert any(n == 2*fib(3*k+2) for k in range(n+1)), (n, cycles)
            print(f'Primitive periodic avoiders: length={n}, words={len(cycles)}, representative={min(cycles)}')
    assert spectrum == [2, 10, 42], spectrum
    print('Fidelity: all 48 Table 1 counts reproduced exactly; a_48=7966')
    print('Complete primitive-period spectrum for 1 <= period <= 48:', spectrum)

    for k in range(1, 5):
        n = 2*fib(3*k+2)
        w = rote_prefix(n)
        assert primitive_root(w) == w
        assert periodic_witness(w) is None, (k, periodic_witness(w))
        print(f'Fidelity: DMRSS Theorem 12, k={k}, period={n}: primitive and avoids')

    a, b = 'S', 'L'
    for i in range(3):
        suitable(a, b)
        aa, bb = len(h(a)), len(h(b))
        assert aa == fib(6*i+5) and bb == fib(6*i+6)
        print(f'Suitable pair level={i}: binary weights=({aa},{bb}), base periods=({2*aa},{2*(aa+2*bb)})')
        if i < 2:
            # Check both terminal families, including rotations of theta(Phi(L)).
            for u in [a, a+b+b]:
                w = h(u) + complement(h(u))
                assert primitive_root(w) == w
                assert periodic_witness(w) is None
        a, b = theta(a), theta(b)

    check_cyclic_descent()

    # The unrelated additional n=0 conjugacy assertion fails, as warned.
    assert rote_prefix(2) == '00' and periodic_witness('01') is None
    print('Scope: additional n=0 prefix-conjugacy assertion fails (00 versus 01); not the target')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
