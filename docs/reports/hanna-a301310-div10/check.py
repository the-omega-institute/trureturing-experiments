#!/usr/bin/env python3
"""Exact finite checks of Hanna's three published divisibility conjectures.

Source: Paul D. Hanna, OEIS A301310/A301311/A301432, NAME, FORMULA,
COMMENTS, and the attributed b-files in this directory. See README.md.
Only integer arithmetic is used; no infinite sum is truncated for evidence.
The finite transforms use the standard Newton/binomial generating identity.
This program is finite evidence, not a proof of a universal assertion.
"""
from hashlib import sha256
from math import comb
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent
LIMIT = 200
MEMBERS = (("301310", 2, 10, False), ("301311", 2, 10, True),
           ("301432", 3, 21, False))


def source_value(k, n, rising):
    if n == 0:
        return 1
    return comb(k*k + n - 1 if rising else k*k, n)


def newton_value(values, m):
    """sum_j Delta^j P(0) m^j, from the Newton binomial expansion."""
    row = values[:]
    result, power = 0, 1
    for _ in range(len(row)):
        result += row[0] * power
        power *= m
        row = [y-x for x, y in zip(row, row[1:])]
    return result


def numerator(values, d):
    """Coefficients of (1-z)^(d+1) sum_{k>=0} P(k) z^k.

    Two coefficients past d are included to test the predicted vanishing.
    This convolution is computed independently of the Newton differences.
    """
    signed = [(-1)**i * comb(d+1, i) for i in range(d+2)]
    return [sum(signed[i] * values[r-i]
                for i in range(min(r, d+1)+1)) for r in range(d+3)]


def divide_z_zplus1(num, d):
    assert num[0] == 0 and num[d+1:] == [0, 0]
    assert sum((-1)**j * c for j, c in enumerate(num)) == 0
    q, previous = [], 0
    for j in range(d-1):
        previous = num[j+1] - previous
        q.append(previous)
    reconstructed = [0] * (d+3)
    for j, c in enumerate(q):
        reconstructed[j+1] += c
        reconstructed[j+2] += c
    assert reconstructed == num
    # Rational reciprocity predicts this exact integer coefficient symmetry.
    padded = num[:d+1] + [0]
    assert padded == list(reversed(padded))
    return q


def homogeneous_eval(coeffs, x, y, degree):
    assert len(coeffs) <= degree+1
    return sum(c * x**j * y**(degree-j) for j, c in enumerate(coeffs))


def read_bfile(aid):
    table = {}
    for line in (ROOT / ("b"+aid+".txt")).read_text().splitlines():
        s = line.strip()
        if s and not s.startswith("#"):
            n, a = map(int, s.split())
            table[n] = a
    assert sorted(table) == list(range(101))
    return table


def main():
    for aid, m, divisor, rising in MEMBERS:
        table = read_bfile(aid)
        values_out, small_q = [1], []
        assert table[0] == 1
        for n in range(1, LIMIT+1):
            d = 2*n
            p = [source_value(k, n, rising) for k in range(d+3)]
            # Newton differences and direct numerator convolution are separate
            # exact evaluations of the polynomial/geometric sum identity.
            a = newton_value(p[:d+1], m)
            num = numerator(p, d)
            assert homogeneous_eval(num[:d+1], m, m+1, d) == a
            q = divide_z_zplus1(num, d)
            quotient = homogeneous_eval(q, m, m+1, d-2)
            assert a == m*(2*m+1)*quotient
            assert m*(2*m+1) == divisor
            assert a % divisor == 0
            if n in table:
                assert a == table[n], (aid, n, a, table[n])
            values_out.append(a)
            if n <= 3:
                small_q.append(q)
        digest = sha256(json.dumps(values_out, separators=(",", ":")).encode()).hexdigest()
        print(f"A{aid}: b-file fidelity n=0..100: 101/101 exact matches")
        print(f"A{aid}: n=1..{LIMIT}: {LIMIT}/{LIMIT} divisible by {divisor}")
        print(f"A{aid}: n=1..{LIMIT}: independent Newton/numerator equality, "
              "tail vanishing, reciprocity, and z(z+1) factorization all checked")
        print(f"A{aid}: a(0..3)={values_out[:4]}; Q(1..3)={small_q}")
        print(f"A{aid}: a({LIMIT}) decimal_digits={len(str(values_out[-1]))}; "
              f"a(0..{LIMIT}) JSON sha256={digest}")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
