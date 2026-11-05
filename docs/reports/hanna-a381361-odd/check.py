#!/usr/bin/env python3
"""Exact finite checks of OEIS A381361's literal defining equation and oddness.

Source: Paul D. Hanna, https://oeis.org/A381361, revision 25.
These finite checks are not a proof of the all-index assertion.
"""
from fractions import Fraction as Q
from math import comb, factorial
from pathlib import Path
import hashlib


def integer_coefficients(last):
    """Triangular Hurwitz multiplication/inversion; includes degrees 0..last."""
    a, c, b = [1], [], []
    for n in range(last + 1):
        c.append(sum(comb(n, k) * a[k] * a[n-k] for k in range(n+1)))
        b.append(1 if n == 0 else
                 -sum(comb(n, k) * c[k] * b[n-k] for k in range(1, n+1)))
        if n < last:
            a.append(sum(comb(n, k) * a[k] * abs(b[n-k])
                         for k in range(n+1)))
    return a, c, b


def ordinary_mul(f, g):
    """Ordinary coefficient convolution, distinct from Hurwitz convolution."""
    assert len(f) == len(g)
    return [sum((f[k] * g[n-k] for k in range(n+1)), Q(0))
            for n in range(len(f))]


def ordinary_inverse(f):
    assert f[0] == 1
    out = [Q(1)]
    for n in range(1, len(f)):
        out.append(-sum((f[k] * out[n-k] for k in range(1, n+1)), Q(0)))
    return out


def zero_primitive(f):
    return [Q(0)] + [f[n-1] / n for n in range(1, len(f))]


def ordinary_exp(f):
    """Taylor exponential by sum f^j/j!, with exactly zero constant f[0]."""
    assert f[0] == 0
    out = [Q(0)] * len(f)
    power = [Q(1)] + [Q(0)] * (len(f)-1)
    for j in range(len(f)):
        out = [out[n] + power[n] / factorial(j) for n in range(len(f))]
        power = ordinary_mul(power, f)
    return out


def literal_map(f):
    reciprocal_square = ordinary_inverse(ordinary_mul(f, f))
    return ordinary_exp(zero_primitive([abs(x) for x in reciprocal_square]))


def main():
    source_path = Path(__file__).with_name("source_terms.tsv")
    source = {}
    for line in source_path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        n, value = map(int, line.split())
        assert n not in source
        source[n] = value
    last = max(source)
    assert set(source) == set(range(last+1))
    assert last == 500
    a, c, b = integer_coefficients(last)
    assert all(a[n] == source[n] for n in range(last+1))
    assert a[:11] == [1, 1, 3, 7, 41, 225, 1435, 11815, 108945,
                       1062145, 12367475]
    assert b[:11] == [1, -2, 0, 16, -64, 0, 2560, -20480, 0,
                       1884160, -21299200]
    assert c[:8] == [1, 2, 8, 32, 192, 1280, 10240, 92160]
    assert all(type(x) is int for row in (a,c,b) for x in row)
    assert all(x % 2 == 1 for x in a)
    assert all(x % 2 == 0 for x in c[1:])
    assert all(x % 2 == 0 for x in b[1:])
    # Check full product identities, including constant coefficients.
    for n in range(last+1):
        assert sum(comb(n,k)*c[k]*b[n-k] for k in range(n+1)) == (n == 0)
    for n in range(last):
        assert a[n+1] == sum(comb(n,k)*a[k]*abs(b[n-k])
                              for k in range(n+1))
        assert (a[n+1]-a[n]) % 2 == 0
    print(f"Source b-file sha256: {hashlib.sha256(source_path.read_bytes()).hexdigest()}")
    print(f"Fidelity: all {last+1} source a(n) agree at n=0..{last}")
    print("Fidelity: displayed source reciprocal-square coefficients agree at n=0..10")
    print(f"Parity: a(0)..a({last}) odd; c(1)..c({last}) and b(1)..b({last}) even")
    print(f"Hurwitz product C*B=1 checked through degree {last}")
    print(f"Derivative recurrence checked at n=0..{last-1}")
    print("a prefix:", a[:11])
    print("c prefix:", c[:8])
    print("b prefix:", b[:11])
    # Independent representation: rational ordinary coefficients, literal map.
    degree = 18
    f = [Q(1)] + [Q(0)] * degree
    for step in range(1, degree+1):
        g = literal_map(f)
        assert g[:step] == f[:step]
        f = g
    assert literal_map(f) == f
    assert [factorial(n)*f[n] for n in range(degree+1)] == a[:degree+1]
    invsq = ordinary_inverse(ordinary_mul(f,f))
    assert [factorial(n)*invsq[n] for n in range(degree+1)] == b[:degree+1]
    primitive = zero_primitive([abs(x) for x in invsq])
    assert primitive[0] == 0
    assert ordinary_exp(primitive) == f
    print(f"Literal rational map exp(integral_0 abs((F^2)^(-1))) stabilizes through degree {degree}")
    print(f"Independent rational coefficients and inverse agree with Hurwitz data at n=0..{degree}")
    print("Scope: finite exact checks only; no all-index proof or Lean probe")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
