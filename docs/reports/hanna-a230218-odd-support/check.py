#!/usr/bin/env python3
"""Independent finite checks for Paul D. Hanna's OEIS A230218 parity conjecture.

Uses Python integers only; the archived OEIS data are checked after reconstruction.
No claim of an all-index proof or Lean verification is made by this program.
"""
from hashlib import sha256
from math import comb
from pathlib import Path

HERE = Path(__file__).resolve().parent
TABLE_SHA256 = '4f6facf9960662a653be9aa37c1a411d64128d405637d20bdfeecb7d100fdc89'


def multiply(f, g, degree):
    out = [0] * (degree + 1)
    for i, v in enumerate(f[:degree + 1]):
        if v:
            for j, w in enumerate(g[:degree + 1 - i]):
                out[i + j] += v * w
    return out


def power(f, exponent, degree):
    out = [1] + [0] * degree
    base = list(f[:degree + 1])
    while exponent:
        if exponent & 1:
            out = multiply(out, base, degree)
        exponent >>= 1
        if exponent:
            base = multiply(base, base, degree)
    return out


def binary_multiply(f, g, mask):
    out = 0
    while g:
        bit = g & -g
        out ^= f << (bit.bit_length() - 1)
        g ^= bit
    return out & mask


def binary_power(f, exponent, degree):
    mask = (1 << (degree + 1)) - 1
    out = 1
    while exponent:
        if exponent & 1:
            out = binary_multiply(out, f, mask)
        exponent >>= 1
        if exponent:
            f = binary_multiply(f, f, mask)
    return out


def main():
    raw = (HERE / 'b230218.txt').read_bytes()
    assert sha256(raw).hexdigest() == TABLE_SHA256
    table = {}
    for line in raw.decode().splitlines():
        if line.strip() and not line.startswith('#'):
            n, a = map(int, line.split())
            assert n not in table
            table[n] = a
    assert set(table) == set(range(521))
    print('OEIS A230218: 521 exact terms, indices 0..520')
    print('OEIS b-file sha256:', TABLE_SHA256)

    # Reconstruct independently with the integer derivative normalization.
    # R=[x^n]F^e, D=[x^(n-1)]F^(e-1)Fprime, N=R-(n+1)D.
    a = [1, 1]
    for n in range(2, 65):
        e = n*n + n + 1
        r = power(a, e, n)[n]
        derivative = [k*a[k] for k in range(1, len(a))]
        d = multiply(power(a, e - 1, n - 1), derivative, n - 1)[n - 1]
        normal = r - (n + 1)*d
        assert n*r == e*d
        assert e*normal == r
        a.append(-normal)
        assert power(a, e, n)[n] == 0
    assert a == [table[n] for n in range(65)]
    assert a[:9] == [1, 1, -3, 14, -85, 504, -4424, 6796, -878157]
    print('Fidelity: reconstructed a(0)..a(64) match the table exactly')
    print('Initial coefficients:', a[:9])
    print('Defining integer rows: zero for n=2..64')

    # A second reconstruction, directly in F_2[[x]], through the full table range.
    bits = 3
    for n in range(2, 521):
        e = n*n + n + 1
        bit = (binary_power(bits, e, n) >> n) & 1
        bits |= bit << n
        assert ((binary_power(bits, e, n) >> n) & 1) == 0
    for n in range(521):
        assert ((bits >> n) & 1) == table[n] % 2
        if n:
            assert table[n] % 2 == int(n & (n - 1) == 0)
    print('Binary reconstruction and OEIS parity agree exactly for n=0..520')
    print('Positive odd support through 520:', [n for n in range(1, 521) if table[n] % 2])

    # Signed Catalan comparison T=x-T^2, C=1+T; no table input.
    t = [0, 1]
    for n in range(2, 41):
        t.append(-sum(t[i]*t[n-i] for i in range(1, n)))
    c = [1] + t[1:]
    for n in range(2, 41):
        e = n*n + n + 1
        row = power(c, e, n)[n]
        rhs = comb(n*n, n) + 2*comb(n*n, n - 1)
        assert row == rhs
        assert row % 2 == 0
    print('Comparison identity: [x^n]C^(n^2+n+1) = binom(n^2,n) + 2 binom(n^2,n-1), n=2..40')
    for n in range(2, 301):
        assert comb(n*n, n) % 2 == 0
    print('Exact binomial parity: even for n=2..300')

    for n in range(2, 100001):
        assert n & (n*n) != n
        if n % 2 == 0:
            p = n & -n  # 2^v2(n)
            assert (n*n) % (p*p) == 0 and n % (p*p) != 0
        else:
            p = (n - 1) & -(n - 1)  # 2^v2(n-1)
            assert n % (2*p) == 1 + p
            assert (n*n) % (2*p) == 1
            assert n % (2*p) not in (0, 1)
    print('Binary obstruction, including both valuation cases: n=2..100000')
    print('Finite exact checks support the proposed general argument; they do not prove it.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
