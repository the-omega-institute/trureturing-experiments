#!/usr/bin/env python3
"""Exact finite checks for OEIS A193237 and A193236 (Bala's two clauses each).

Source: Paul D. Hanna's OEIS definitions (2011); Peter Bala's %F conjectures
(2021), revisions 18 and 20 respectively. Prefixes below quote %S/%T/%U.
OEIS data are subject to https://oeis.org/LICENSE.
No finite check here establishes the unbounded conjectures.
"""
import argparse
from fractions import Fraction
import hashlib
from pathlib import Path

LIMIT = 3505
PREFIXES = {
    'A193237': [1,3,18,141,1125,9261,78255,673137,5864238,51592770,
                457484382,4082618376,36627119109,330070717935,2985857903655,
                27099681108948,246666402397287,2250904657271427,
                20586440729350197,188659279149885810,1732045683183434379],
    'A193236': [1,3,18,111,765,5481,40581,306099,2342034,18108270,
                141176412,1108011312,8744143401,69325981191,551800999215,
                4406974587918,35300439813735,283495238613855,
                2281964065354899,18406084773140820,148734744069134439],
}


def source_series(sid, limit):
    f = [0] * (limit + 1)
    j = 0
    while j * (j + 1) // 2 <= limit:
        t = j * (j + 1) // 2
        f[t] = ((2*j + 1) * 3**j * (-1)**t if sid == 'A193237'
                else (-3)**j * (2*j + 1))
        j += 1
    return f


def reciprocal_cube_root(f):
    # Coefficients of 3 F A' + F' A = 0, A(0)=1, over Q.
    # Exact divisibility is checked at every step, never rounded.
    a = [1]
    sparse = [(k, v) for k, v in enumerate(f) if k and v]
    for n in range(1, len(f)):
        numerator = -sum((3*n - 2*k)*v*a[n-k]
                         for k, v in sparse if k <= n)
        q, r = divmod(numerator, 3*n)
        assert r == 0, ('nonintegral coefficient', n, numerator, 3*n)
        a.append(q)
    return a


def mul(a, b, limit, modulus=None):
    out = [0] * (limit + 1)
    aa = [(i, x) for i, x in enumerate(a[:limit+1]) if x]
    bb = [(j, y) for j, y in enumerate(b[:limit+1]) if y]
    for i, x in aa:
        for j, y in bb:
            if i+j > limit:
                break
            out[i+j] += x*y
    if modulus:
        out = [x % modulus for x in out]
    return out


def binomial_construction(f, limit):
    # Independent binomial sum (1+9H)^(-1/3), through degree limit.
    assert f[0] == 1 and all(v % 9 == 0 for v in f[1:])
    h = [0] + [v//9 for v in f[1:limit+1]]
    power = [1] + [0]*limit
    out = [0]*(limit+1)
    c = Fraction(1)
    for k in range(limit+1):
        assert c.denominator == 1, ('nonintegral binomial coefficient', k)
        for n in range(limit+1):
            out[n] += c.numerator*power[n]
        power = mul(power, h, limit)
        c *= Fraction(-3*(3*k+1), k+1)
    return out


def verify_frobenius_transfer(f, a, p, exponent):
    limit = len(a)-1
    ff = [x % p for x in f]
    support = [1] + [0]*limit
    for _ in range(exponent):
        support = mul(support, ff, limit, p)
    b = [x % p for x in a[:limit//p+1]]
    if p == 5:
        b = mul(b, b, limit//p, p)
    expanded = [b[n//p] if n % p == 0 else 0 for n in range(limit+1)]
    rhs = mul(expanded, support, limit, p)
    assert rhs == [x % p for x in a], ('Frobenius transfer', p)
    missing = 4 if p == 5 else 5
    assert all(support[n] == 0 for n in range(missing, limit+1, p))


def check_bfile(directory, sid, a):
    path = directory / ('b' + sid[1:] + '.txt')
    raw = path.read_bytes()
    rows = [tuple(map(int, line.split())) for line in raw.decode().splitlines()
            if line.strip() and not line.lstrip().startswith('#')]
    assert [n for n, _ in rows] == list(range(1001)), (sid, 'b-file indices')
    assert all(a[n] == v for n, v in rows), (sid, 'b-file mismatch')
    print(f'{sid}: b-file indices 0..1000 match exactly; '
          f'sha256={hashlib.sha256(raw).hexdigest()}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bfile-dir', type=Path,
                        help='optional directory containing the two OEIS b-files')
    args = parser.parse_args()
    assert [j*(j+1)//2 % 5 for j in range(5)] == [0,1,3,1,0]
    assert [j*(j+1)//2 % 7 for j in range(7)] == [0,1,3,6,3,1,0]
    print('Triangular residue tables at 5 and 7: exact match')
    for sid, prefix in PREFIXES.items():
        f = source_series(sid, LIMIT)
        a = reciprocal_cube_root(f)
        assert a[:len(prefix)] == prefix, (sid, 'OEIS prefix mismatch')
        print(f'{sid}: source prefix indices 0..20 reproduced exactly')
        assert a[:101] == binomial_construction(f, 100)
        unit = mul(mul(mul(a, a, 100), a, 100), f, 100)
        assert unit == [1]+[0]*100, (sid, 'A^3 F equation')
        print(f'{sid}: independent binomial construction and A^3 F=1 through degree 100')
        if args.bfile_dir:
            check_bfile(args.bfile_dir, sid, a)
        for p, residue, exponent in [(5,4,3),(7,5,2)]:
            indices = list(range(residue, LIMIT+1, p))
            assert all(a[m] % p == 0 for m in indices), (sid, p)
            # n=0..200 includes index 1004 (p=5) and 1405 (p=7).
            assert all(a[p*n+residue] % p == 0 for n in range(201))
            verify_frobenius_transfer(f, a, p, exponent)
            print(f'{sid}: mod {p}, residue {residue}: {len(indices)} cases, '
                  f'n=0..{len(indices)-1}, indices {indices[0]}..{indices[-1]}; '
                  'no violations; full transfer identity through degree 3505')
        print(f'{sid}: a(4)={a[4]}=5*{a[4]//5}; a(5)={a[5]}=7*{a[5]//7}')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
