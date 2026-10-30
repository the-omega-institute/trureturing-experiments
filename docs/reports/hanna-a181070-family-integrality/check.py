#!/usr/bin/env python3
"""Exact finite checks for the three Hanna variable-power exponentials.

This program tests finite ranges; it does not prove the conjectures.
The adjacent OEIS b-files are fidelity fixtures, not integrality evidence.
"""

from fractions import Fraction
from hashlib import sha256
from math import comb, factorial, isqrt
from pathlib import Path


DEGREE = 140
DIRECT_DEGREE = 18
PRIMES = tuple(p for p in range(2, DEGREE + 1)
               if all(p % d for d in range(2, isqrt(p) + 1)))
PREFIXES = {
    "A181070": (1, 1, 2, 4, 8, 23, 88, 379, 3044, 32116),
    "A181074": (1, 1, 2, 5, 23, 231, 5405, 322799, 42761356, 12597156231),
    "A181080": (1, 1, 2, 4, 14, 83, 774, 10641, 255918, 14643874),
}
LOG_PREFIXES = {
    "A181070": (0, 1, 3, 7, 15, 66, 357, 1891),
    "A181074": (0, 1, 3, 10, 71, 1026, 30912),
    "A181080": (0, 1, 3, 7, 39, 336, 4077),
}


def vp(p, n):
    assert n > 0
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def logarithm(seq, degree, missing_extra_factor=False):
    """Literal source formula, coefficientwise over Q."""
    out = [Fraction(0) for _ in range(degree + 1)]
    for n in range(1, degree + 1):
        last_k = degree - n if seq == "A181074" else min(n, degree - n)
        for k in range(last_k + 1):
            b = comb(n + k - 1, k) if seq == "A181074" else comb(n, k)
            e = n - k + 1 if seq == "A181080" else k + 1
            e -= int(missing_extra_factor)
            out[n + k] += Fraction(b ** e, n)
    return out


def exp_via_recurrence(log):
    """a_0=1 and m*a_m=sum_{j=1}^m j*L_j*a_{m-j}."""
    a = [Fraction(1)]
    for m in range(1, len(log)):
        a.append(sum((j * log[j] * a[m - j] for j in range(1, m + 1)),
                     Fraction(0)) / m)
    return a


def mul(a, b, degree):
    out = [Fraction(0) for _ in range(degree + 1)]
    for i in range(min(len(a), degree + 1)):
        for j in range(min(len(b), degree - i + 1)):
            out[i + j] += a[i] * b[j]
    return out


def exp_via_taylor(log, degree):
    """Independent truncated sum of L^r/r!, with L_0=0."""
    assert log[0] == 0
    power = [Fraction(1)] + [Fraction(0)] * degree
    out = [Fraction(0)] * (degree + 1)
    for r in range(degree + 1):
        for m in range(degree + 1):
            out[m] += power[m] / factorial(r)
        power = mul(power, log[:degree + 1], degree)
    return out


def check_entry(seq):
    log = logarithm(seq, DEGREE)
    c = [m * log[m] for m in range(DEGREE + 1)]
    assert tuple(c[:len(LOG_PREFIXES[seq])]) == LOG_PREFIXES[seq], seq
    assert all(v.denominator == 1 for v in c), seq
    a = exp_via_recurrence(log)
    assert tuple(a[:10]) == PREFIXES[seq], seq
    assert all(v.denominator == 1 for v in a), seq
    assert a[:DIRECT_DEGREE + 1] == exp_via_taylor(log, DIRECT_DEGREE), seq
    fixture = Path(__file__).with_name("b" + seq[1:] + ".txt")
    rows = [line.split() for line in fixture.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")]
    for i, text in rows:
        assert a[int(i)] == int(text), (seq, i)
    dold_count = 0
    for p in PRIMES:
        for m in range(1, DEGREE // p + 1):
            modulus = p ** (1 + vp(p, m))
            assert (c[p * m] - c[m]) % modulus == 0, (seq, p, m)
            dold_count += 1
    # The source's comparison formula omits the extra +1 in the exponent.
    comparison = exp_via_recurrence(logarithm(seq, 20, True))
    obstruction = next((m, v) for m, v in enumerate(comparison)
                       if v.denominator != 1)
    print(f"{seq}: first10={','.join(str(v) for v in a[:10])}")
    print(f"{seq}: fidelity b-file indices={rows[0][0]}..{rows[-1][0]}, "
          f"rows={len(rows)}, sha256={sha256(fixture.read_bytes()).hexdigest()}")
    print(f"{seq}: integral coefficients and integral m*[x^m]L for m=0..{DEGREE}; "
          f"Taylor/recurrence agreement m=0..{DIRECT_DEGREE}; "
          f"Dold pairs={dold_count} with p*m<={DEGREE}")
    print(f"{seq}: coefficient[{DEGREE}] sha256="
          f"{sha256(str(a[-1]).encode()).hexdigest()}")
    print(f"{seq}: comparison without +1 first noninteger="
          f"coefficient[{obstruction[0]}]={obstruction[1]}")
    if seq == "A181080":
        # Reindex by j=n-k; retain the original source series above.
        reindexed = [Fraction(0) for _ in range(DEGREE + 1)]
        for n in range(1, DEGREE + 1):
            for j in range(n + 1):
                m = 2 * n - j
                if m <= DEGREE:
                    reindexed[m] += Fraction(comb(n, j) ** (j + 1), n)
        assert reindexed == log
        print("A181080: j=n-k, weight 2*n-j reindexing agrees through degree 140")


def check_variable_power_congruence():
    primes = [p for p in PRIMES if p <= 19]
    cases = units = multiples = zeros = 0
    for family in ("binomial", "negative-binomial"):
        for n in range(1, 41):
            for k in range(41):
                if family == "binomial":
                    b = comb(n, k) if k <= n else 0
                else:
                    b = comb(n + k - 1, k)
                assert k * b % n == 0, (family, n, k)
                for p in primes:
                    q = p ** (1 + vp(p, n))
                    if family == "binomial":
                        bp = comb(p * n, p * k) if k <= n else 0
                    else:
                        bp = comb(p * (n + k) - 1, p * k)
                    assert (bp - b) % q == 0, (family, n, k, p, "base")
                    assert (pow(bp, p * k + 1, q) - pow(b, k + 1, q)) % q == 0, (
                        family, n, k, p, "variable-power")
                    cases += 1
                    if b == 0:
                        zeros += 1
                    elif b % p:
                        units += 1
                    else:
                        multiples += 1
    print("Variable-power congruences: families=binomial,negative-binomial; "
          "n=1..40; k=0..40; primes=2,3,5,7,11,13,17,19; "
          f"cases={cases}; units={units}; nonzero-p-multiples={multiples}; zeros={zeros}")


def main():
    for seq in PREFIXES:
        check_entry(seq)
    check_variable_power_congruence()
    print("Finite exact checks only; the all-degree statements are unproved by this program.")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
