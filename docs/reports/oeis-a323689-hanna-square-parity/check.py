#!/usr/bin/env python3
"""Exact finite fidelity and parity checks for Paul D. Hanna's two OEIS series.

All Gaussian integers are pairs of Python integers. No floating-point arithmetic
or external algebra package is used. Run from any directory with Python 3.
"""
from hashlib import sha256
from math import comb, isqrt
from pathlib import Path

HERE = Path(__file__).resolve().parent
ZERO = (0, 0)
ONE = (1, 0)
IPOW = (ONE, (0, 1), (-1, 0), (0, -1))


def add(a, b):
    return (a[0] + b[0], a[1] + b[1])


def sub(a, b):
    return (a[0] - b[0], a[1] - b[1])


def mul(a, b):
    return (a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0])


def scale(n, a):
    return (n*a[0], n*a[1])


def ipow(n):
    return IPOW[n % 4]


def rational_series_division(m, degree):
    """Divide (t+i)^m by (1+i*t)^(m+1), using its unit constant term."""
    numerator = [scale(comb(m, j), ipow(m-j)) for j in range(m+1)]
    denominator = [scale(comb(m+1, j), ipow(j)) for j in range(m+2)]
    result = []
    for k in range(degree+1):
        value = numerator[k] if k <= m else ZERO
        for j in range(1, min(m+1, k)+1):
            value = sub(value, mul(denominator[j], result[k-j]))
        result.append(value)
    return result


def direct_coefficient(m, k):
    """Binomial expansion of numerator and reciprocal, without Delannoy numbers."""
    value = ZERO
    for j in range(min(m, k)+1):
        scalar = comb(m, j)*comb(m+k-j, k-j)
        value = add(value, scale(scalar, ipow(m-j-(k-j))))
    return value


def delannoy(m, k):
    return sum((1 << j)*comb(m, j)*comb(k, j) for j in range(min(m, k)+1))


def delannoy_coefficient(m, k):
    return scale(delannoy(m, k), ipow(m-k))


def divisors(n):
    for d in range(1, isqrt(n)+1):
        if n % d == 0:
            yield d
            if d*d != n:
                yield n//d


def coefficient(anum, n, kernel):
    if anum == 'A323689':
        if n == 0:
            return ONE
        return sum_pairs(kernel(d, n//d) for d in divisors(n))
    return sum_pairs(kernel(d-1, (n+1)//d-1) for d in divisors(n+1))


def sum_pairs(values):
    value = ZERO
    for term in values:
        value = add(value, term)
    return value


def literal_truncation(anum, degree):
    """Truncate the source's actual rational summands by polynomial division."""
    out = [ZERO]*(degree+1)
    if anum == 'A323689':
        out[0] = ONE
        for m in range(1, degree+1):
            q = rational_series_division(m, degree//m)
            for k in range(1, len(q)):
                out[m*k] = add(out[m*k], q[k])
    else:
        for m in range(degree+1):
            q = rational_series_division(m, (degree-m)//(m+1))
            for k, v in enumerate(q):
                n = m+(m+1)*k
                out[n] = add(out[n], v)
    return out


def read_table(anum):
    path = HERE / ('b'+anum[1:]+'.txt')
    rows = {}
    for line in path.read_text().splitlines():
        line = line.split('#', 1)[0].strip()
        if line:
            n, a = map(int, line.split())
            assert n not in rows
            rows[n] = a
    assert sorted(rows) == list(range(len(rows)))
    print(f'{anum} source table: rows={len(rows)} range=0..{len(rows)-1} '
          f'sha256={sha256(path.read_bytes()).hexdigest()}')
    return rows


def main():
    # A third, independent Delannoy computation by the lattice-path recurrence.
    size = 24
    grid = [[1]*(size+1) for _ in range(size+1)]
    for m in range(1, size+1):
        for k in range(1, size+1):
            grid[m][k] = grid[m-1][k]+grid[m][k-1]+grid[m-1][k-1]
    for m in range(size+1):
        divided = rational_series_division(m, size)
        for k in range(size+1):
            d = delannoy(m, k)
            assert grid[m][k] == d == delannoy(k, m)
            assert d % 2 == 1
            assert divided[k] == direct_coefficient(m, k) == delannoy_coefficient(m, k)
            pair = add(delannoy_coefficient(m, k), delannoy_coefficient(k, m))
            assert pair[1] == 0 and pair[0] % 2 == 0
    print('Kernel check: m,k=0..24; 625 exact identities; division, binomial '
          'expansion and Delannoy recurrence agree')

    for anum, maximum in [('A323689', 10000), ('A324300', 10201)]:
        table = read_table(anum)
        assert max(table) == maximum
        literal = literal_truncation(anum, 80)
        for n in range(81):
            assert literal[n] == coefficient(anum, n, direct_coefficient)
            assert literal[n] == (table[n], 0)
        print(f'{anum} fidelity: literal rational summands by polynomial '
              'division reproduce published a(0)..a(80)')
        for n, expected in table.items():
            value = coefficient(anum, n, direct_coefficient)
            via_delannoy = coefficient(anum, n, delannoy_coefficient)
            assert value == via_delannoy == (expected, 0), (anum, n, value, expected)
            target = n if anum == 'A323689' else n+1
            assert (value[0] % 2 == 1) == (isqrt(target)**2 == target), (anum, n)
        print(f'{anum} full-table fidelity and expected parity: '
              f'{len(table)} cases; coefficient, integrality and parity mismatches=0')
        print(f'{anum} first 17 coefficients: {[table[n] for n in range(17)]}')
    assert delannoy(4, 4) == 321 and delannoy(2, 8) == 145
    assert coefficient('A323689', 16, direct_coefficient) == (321-2*145, 0)
    assert coefficient('A324300', 8, direct_coefficient) == (2+13, 0)
    print('Examples: A323689(16)=321-2*145=31; A324300(8)=2+13=15')
    print('Finite checks support the transcription and proposed route; they do not prove either unbounded claim.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
