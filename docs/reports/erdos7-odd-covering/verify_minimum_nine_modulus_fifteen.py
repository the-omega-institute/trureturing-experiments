#!/usr/bin/env python3
"""Verify the absence of an odd distinct cover omitting 3, 5, 7 and 15."""

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
    count, total = map(int, result.stdout.split())
    if count != sum(primes[11:PREFIX + 1]):
        raise ValueError("prefix prime count disagrees with the continuation sieve")
    if not 500 * total < 499 * PREFIX_SCALE:
        raise ValueError("prefix reserve is too small")
    print(f"Tail threshold {upper_tail}/{SCALE}, prefix seed {seed}/{PREFIX_SCALE}")
    print(f"{count} stages at primes 11 through {PREFIX}")
    print(f"Total {total}/{PREFIX_SCALE}, reserve above 1/500")


if __name__ == "__main__":
    main()
