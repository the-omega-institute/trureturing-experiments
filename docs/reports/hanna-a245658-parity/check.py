#!/usr/bin/env python3
"""Exact finite checks for Hanna's OEIS A245658 parity conjecture.

Source: A245658 revision 24 and A247031 revision 12, both 2014.
The source's binomial-cube triangular transform is evaluated directly.
No finite check here proves the all-index conjecture.
"""
from math import comb

LIMIT = 100


def multiply(p, q, limit):
    out = [0] * (limit + 1)
    for i, pi in enumerate(p[:limit + 1]):
        for j, qj in enumerate(q[:limit + 1 - i]):
            out[i + j] += pi * qj
    return out


def power(p, exponent, limit):
    result = [1] + [0] * limit
    for _ in range(exponent):
        result = multiply(result, p, limit)
    return result


def kernel(r, limit):
    return [(-1)**j * comb(r, j)**3 for j in range(min(r, limit) + 1)]


def main():
    nmax = LIMIT
    kernels = [power(kernel(m + 1, nmax - m), 3, nmax - m)
               for m in range(nmax + 1)]
    b = [1]
    for n in range(1, nmax + 1):
        b.append(-sum(b[m] * kernels[m][n - m] for m in range(n)))
    a = [1]
    for n in range(1, nmax + 1):
        lower = power(a, 3, n)[n]
        quotient, remainder = divmod(b[n] - lower, 3)
        assert remainder == 0, ('nonintegral cube root', n, remainder)
        a.append(quotient)
    expected_a = [1, 1, 22, 1624, 264962, 79136637, 38671111558,
                  28642761340956, 30413158977739302, 44371247722932580948,
                  86030010617294195269924, 215785986252313362542154058,
                  684854225048414942925120331598,
                  2699657503162253569920254747627596,
                  13010506207186236974375590663943378970]
    expected_b = [1, 3, 69, 5005, 806148, 239220375, 116532061510,
                  86173621173099, 91417549409916684, 133300597778263476112,
                  258360728839130761571757, 647880493609691058921741273,
                  2055869510173976408422116133220,
                  8103111707775918586405906798540650,
                  39047811321420953231675462397758519802]
    assert a[:15] == expected_a
    assert b[:15] == expected_b
    print('FIDELITY: both OEIS displayed coefficient lists agree at indices 0..14')
    print('A245658 initial coefficients:', a[:7])
    assert power(a, 3, nmax) == b
    residual = [sum(b[m] * kernels[m][n - m] for m in range(n + 1))
                for n in range(nmax + 1)]
    assert residual == [1] + [0] * nmax
    print('SOURCE: triangular transform = 1 and A^3 = B through degree', nmax)
    c = []
    for n in range(nmax + 1):
        cn, rem = divmod(comb(4*n, n), 3*n + 1)
        assert rem == 0
        c.append(cn)
    c4 = power(c, 4, nmax)
    assert c == [1] + c4[:nmax]
    c3 = power(c, 3, nmax)
    assert all((b[n] - c3[n]) % 9 == 0 for n in range(nmax + 1))
    for m, p3 in enumerate(kernels):
        r = m + 1
        assert all((v - (-1)**j * comb(3*r, j)) % 9 == 0
                   for j, v in enumerate(p3))
    print('INTEGRALITY ROUTE: kernel congruence and B = C^3 modulo 9 verified through degree', nmax)
    y = [v % 2 for v in a]
    y4 = power(y, 4, nmax)
    assert all(y[n] == (1 if n == 0 else y4[n - 1] % 2)
               for n in range(nmax + 1))
    assert all(y4[n] % 2 == (y[n//4] if n % 4 == 0 else 0)
               for n in range(nmax + 1))
    expected_odd = []
    four_power = 1
    while (four_power - 1)//3 <= nmax:
        expected_odd.append((four_power - 1)//3)
        four_power *= 4
    odd = [n for n, v in enumerate(a) if v % 2 == 1]
    assert odd == expected_odd
    print('PARITY: Y = 1 + xY^4 and fourth-power Frobenius agree through degree', nmax)
    print('ODD INDICES:', odd)
    print('Expected settlement remains an all-index proof obligation; this is finite evidence.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
