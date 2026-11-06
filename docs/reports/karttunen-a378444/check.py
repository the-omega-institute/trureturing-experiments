#!/usr/bin/env python3
"""Exact finite checks for Karttunen's OEIS A378444 ceiling conjecture.

Uses only the Python standard library. The theorem is not proved by this check.
"""
import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import product
from math import gcd, prod
from pathlib import Path
from urllib.request import Request, urlopen

LIMIT = 65537
BFILES = {
    "A378444": "5df9ecbacda6bd7ee7ee6af0c49056b8cea3b99c1667603ad1d42700aa420939",
    "A174273": "1a29482ee696937cf440f64ad9c950ce2ccd3cbf5898a605bffbe05721b4c88a",
}


def bfile(sequence, directory):
    filename = "b" + sequence[1:] + ".txt"
    url = "https://oeis.org/" + sequence + "/" + filename
    if directory is None:
        with urlopen(Request(url, headers={"User-Agent": "Mozilla/5.0"}), timeout=60) as response:
            data = response.read()
    else:
        data = (directory / filename).read_bytes()
    digest = sha256(data).hexdigest()
    assert digest == BFILES[sequence], (sequence, "source bytes changed", digest)
    rows = {}
    for line in data.decode("utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        n, value = map(int, line.split())
        assert n not in rows, (sequence, "duplicate index", n)
        rows[n] = value
    assert sorted(rows) == list(range(1, LIMIT + 1)), (sequence, "unexpected range")
    print("SOURCE", sequence, "rows=" + str(len(rows)), "sha256=" + digest)
    return rows


def factorizations(limit):
    smallest = list(range(limit + 1))
    for p in range(2, limit + 1):
        if smallest[p] == p:
            for m in range(p * p, limit + 1, p):
                if smallest[m] == m:
                    smallest[m] = p
    factors = [()] * (limit + 1)
    for n in range(2, limit + 1):
        m, row = n, []
        while m > 1:
            p, e = smallest[m], 0
            while m % p == 0:
                m //= p
                e += 1
            row.append((p, e))
        factors[n] = tuple(row)
    return factors


def ratio(factors):
    return sum((Fraction(e, p) for p, e in factors), Fraction(0, 1))


def check_large_factor_boxes():
    cases = [(), ((2, 1000),), ((3, 120), (5, 80)),
             ((2, 17), (3, 8), (5, 6)), ((2, 20), (3, 9), (7, 4)),
             ((2, 12), (3, 4), (5, 2), (7, 6)), ((2, 2), (65537, 8))]
    for factors in cases:
        c = b = 0
        for exponents in product(*(range(e + 1) for _, e in factors)):
            r = ratio(tuple((p, e) for (p, _), e in zip(factors, exponents)))
            c += int(r.numerator % 2 == 0)
            a = next((e for (p, _), e in zip(factors, exponents) if p == 2), 0)
            b += int(a % 2 == 0)
        alpha = next((e for p, e in factors if p == 2), 0)
        expected_b = (alpha // 2 + 1) * prod(e + 1 for p, e in factors if p != 2)
        assert b == expected_b and c == (b + 1) // 2, (factors, c, b)
    print("LARGE_FACTOR_BOXES", "checked=" + str(len(cases)), "mismatches=0")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bfile-dir", type=Path,
                        help="optional directory of the two pinned OEIS b-files")
    args = parser.parse_args()
    published = {name: bfile(name, args.bfile_dir) for name in BFILES}
    factors = factorizations(LIMIT)
    numerator_even = [False] * (LIMIT + 1)
    two_exponent_even = [False] * (LIMIT + 1)
    for d in range(1, LIMIT + 1):
        r = ratio(factors[d])
        derivative = sum(e * (d // p) for p, e in factors[d])
        assert r == Fraction(derivative, d), ("logarithmic derivative", d)
        assert r.numerator == derivative // gcd(d, derivative), ("reduction", d)
        numerator_even[d] = r.numerator % 2 == 0
        alpha = next((e for p, e in factors[d] if p == 2), 0)
        odd_omega = sum(e for p, e in factors[d] if p != 2)
        two_exponent_even[d] = alpha % 2 == 0
        predicted = alpha % 2 == 0 and (alpha // 2 + odd_omega) % 2 == 0
        assert numerator_even[d] == predicted, ("numerator parity", d, r)
    print("EXACT_RATIONAL_CHECKS", "range=1..65537", "mismatches=0")

    # Direct divisor sums from exact rational numerators; no conjectural formula here.
    c, b = [0] * (LIMIT + 1), [0] * (LIMIT + 1)
    for d in range(1, LIMIT + 1):
        for n in range(d, LIMIT + 1, d):
            c[n] += int(numerator_even[d])
            b[n] += int(two_exponent_even[d])
    for n in range(1, LIMIT + 1):
        alpha = next((e for p, e in factors[n] if p == 2), 0)
        box_size = (alpha // 2 + 1) * prod(e + 1 for p, e in factors[n] if p != 2)
        assert b[n] == box_size, ("box size", n, b[n], box_size)
        assert c[n] == (b[n] + 1) // 2, ("conjecture", n, c[n], b[n])
        assert c[n] == published["A378444"][n], ("A378444 fidelity", n)
        assert b[n] == published["A174273"][n], ("A174273 fidelity", n)
    print("BFILE_FIDELITY", "A378444=65537", "A174273=65537", "mismatches=0")
    print("CONJECTURE_AND_BOX_FORMULA", "range=1..65537", "mismatches=0")
    # An ancillary A174273 formula says a(2n+1)=sigma(2n+1). It is not used.
    assert b[3] == 2 and sum(d for d in range(1, 4) if 3 % d == 0) == 4
    print("SOURCE_ANCILLARY_DISCREPANCY", "A174273(3)=2", "sigma(3)=4",
          "definition_and_bfile_used")
    check_large_factor_boxes()
    for n in [1, 2, 4, 8, 9, 12, 16, 36, 60, 64, 144, 65536, 65537]:
        print("EXAMPLE", n, "A378444=" + str(c[n]), "A174273=" + str(b[n]))
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
