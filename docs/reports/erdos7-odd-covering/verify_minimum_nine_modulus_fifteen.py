#!/usr/bin/env python3
"""Verify modulus and residue constraints for a distinct odd cover of minimum nine."""

from fractions import Fraction
from pathlib import Path
import subprocess
import tempfile

from verify_rh_prime_continuation import SCALE, ceiling, prime_sieve, verify_quadratic

PREFIX = 80_021
PREFIX_SCALE = 10**14


def main():
    primes = prime_sieve(5_000_000)
    upper_tail = verify_quadratic(primes, capture=PREFIX)
    seed = ceiling(Fraction(PREFIX_SCALE * upper_tail, SCALE))
    source = Path(__file__).with_suffix(".cpp").resolve()
    with tempfile.TemporaryDirectory(prefix="erdos-7-modulus-fifteen-") as directory:
        binary = Path(directory) / "verify"
        subprocess.run(
            ["c++", "-std=c++17", "-O3", "-Wall", "-Wextra",
             "-fsanitize=undefined", "-fno-sanitize-recover=undefined",
             str(source), "-o", str(binary)],
            check=True,
        )
        result = subprocess.run(
            [str(binary)], input=f"{seed}\n", text=True,
            capture_output=True, check=True,
        )
    count, total, missing_both, suffix_floor = map(int, result.stdout.split())
    if count != sum(primes[11:PREFIX + 1]):
        raise ValueError("prefix prime count disagrees with the continuation sieve")
    if not 500 * total < 499 * PREFIX_SCALE:
        raise ValueError("prefix reserve is too small")
    credit = Fraction(3, 100)
    if Fraction(suffix_floor, PREFIX_SCALE) < credit:
        raise ValueError("suffix floor does not supply the claimed credit")
    other_five_slots = Fraction(7, 95)
    aligned_fifteen = Fraction(4, 19) * Fraction(4, 15)
    fifteen_with_no_active_forty_five = Fraction(4, 19) * (Fraction(2, 5) - Fraction(2, 15))
    bound = (Fraction(missing_both, PREFIX_SCALE) - credit * other_five_slots
             + (1 - credit) * max(aligned_fifteen, fifteen_with_no_active_forty_five))
    if not bound < Fraction(999, 1000):
        raise ValueError("residue constraints have insufficient reserve")
    print(f"Tail threshold {upper_tail}/{SCALE}, prefix seed {seed}/{PREFIX_SCALE}")
    print(f"{count} stages at primes 11 through {PREFIX}")
    print(f"Total {total}/{PREFIX_SCALE}, reserve above 1/500")
    print(f"Missing 15 and 21 comparison {missing_both}/{PREFIX_SCALE}")
    print(f"Suffix floor {suffix_floor}/{PREFIX_SCALE}")
    print(f"Residue constraints without 21 bound {bound}, reserve above 1/1000")


if __name__ == "__main__":
    main()
