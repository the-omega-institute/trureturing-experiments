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
    (count, total, missing_both, missing_three, suffix_floor, zero_tail, limited_ternary_union,
     no_twenty_one_or_square, no_small_powers, early_mixed_bound) = map(int, result.stdout.split())
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
    arbitrary_fifteen = Fraction(8, 95)
    without_twenty_one_and_thirty_three = (
        Fraction(missing_three, PREFIX_SCALE) - credit * other_five_slots
        + (1 - credit) * arbitrary_fifteen
    )
    if not without_twenty_one_and_thirty_three < Fraction(1999, 2000):
        raise ValueError("missing 21 and 33 comparison has insufficient reserve")
    full_five_charge = Fraction(3, 19)
    zero_tail_with_credit = Fraction(zero_tail, PREFIX_SCALE) - credit * full_five_charge
    if not zero_tail_with_credit > Fraction(101, 100):
        raise ValueError("zero-tail comparison floor does not exceed 1.01")
    if not 2000 * limited_ternary_union < 1999 * PREFIX_SCALE:
        raise ValueError("limited ternary union has insufficient reserve")
    if not 50_000 * max(no_twenty_one_or_square, no_small_powers) < 49_999 * PREFIX_SCALE:
        raise ValueError("missing-square comparisons have insufficient reserve")
    charge_at_seven = 1 - Fraction(81, 10_000)
    charge_at_five = 1 - Fraction(3486, 100_000)
    no_twenty_one = Fraction(early_mixed_bound, PREFIX_SCALE) - charge_at_seven * Fraction(12, 205)
    no_twenty_one_or_thirty_five = no_twenty_one - charge_at_seven * Fraction(24, 779)
    no_sixty_three = no_twenty_one - charge_at_seven * Fraction(4, 205)
    no_sixty_three_or_one_hundred_five = no_sixty_three - charge_at_seven * Fraction(48, 3895)
    smaller_fifteen = no_sixty_three - charge_at_five * Fraction(4, 19) * (Fraction(2, 5) - Fraction(9, 25))
    if not max(no_twenty_one_or_thirty_five, no_sixty_three_or_one_hundred_five,
               smaller_fifteen) < Fraction(9995, 10_000):
        raise ValueError(
            f"early mixed losses have insufficient reserve: {early_mixed_bound}, "
            f"{no_twenty_one_or_thirty_five}, {no_sixty_three_or_one_hundred_five}, {smaller_fifteen}"
        )
    overlap_of_fifteen_and_thirty_five = Fraction(48, 3895)
    overlap_of_fifteen_and_sixty_three = Fraction(4, 205) * Fraction(4, 19)
    aligned_without_sixty_three = no_sixty_three - charge_at_seven * overlap_of_fifteen_and_thirty_five
    aligned_without_one_hundred_five = no_twenty_one - charge_at_seven * (
        Fraction(48, 3895) + overlap_of_fifteen_and_thirty_five
        + overlap_of_fifteen_and_sixty_three
    )
    if not max(aligned_without_sixty_three, aligned_without_one_hundred_five) < Fraction(999, 1000):
        raise ValueError("phase overlap constraints have insufficient reserve")
    print(f"Tail threshold {upper_tail}/{SCALE}, prefix seed {seed}/{PREFIX_SCALE}")
    print(f"{count} stages at primes 11 through {PREFIX}")
    print(f"Total {total}/{PREFIX_SCALE}, reserve above 1/500")
    print(f"Missing 15 and 21 comparison {missing_both}/{PREFIX_SCALE}")
    print(f"Suffix floor {suffix_floor}/{PREFIX_SCALE}")
    print(f"Residue constraints without 21 bound {bound}, reserve above 1/1000")
    print(f"Missing 15, 21 and 33 comparison {missing_three}/{PREFIX_SCALE}")
    print(f"Missing 21 and 33 bound {without_twenty_one_and_thirty_three}, reserve above 1/2000")
    print(f"All slots, zero tail and prime-5 credit give comparison at least {zero_tail_with_credit}")
    print(f"Missing 21 and ternary union at most 5/36 give bound {limited_ternary_union}/{PREFIX_SCALE}")
    print(f"Missing 21 and 25 give bound {no_twenty_one_or_square}/{PREFIX_SCALE}")
    print(f"Missing 25, 27 and 49 give bound {no_small_powers}/{PREFIX_SCALE}")
    print(f"All early mixed slots give bound {early_mixed_bound}/{PREFIX_SCALE}")
    print(f"Missing 21 and 35 give bound {no_twenty_one_or_thirty_five}")
    print(f"Missing 21, 63 and 105 give bound {no_sixty_three_or_one_hundred_five}")
    print(f"Missing 21 and 63 with ternary 15 mass at most 9/25 gives bound {smaller_fifteen}")
    print(f"Missing 21 and 63 with 15 and 35 aligned modulo 5 gives bound {aligned_without_sixty_three}")
    print(f"Missing 21 and 105 with both phase alignments gives bound {aligned_without_one_hundred_five}")


if __name__ == "__main__":
    main()
