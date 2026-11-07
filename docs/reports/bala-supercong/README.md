# Bala diagonal Euler-product congruences

Run `python3 docs/reports/bala-supercong/check.py` from the repository root.
The expected final line is `ALL_CHECKS_OK`; the checked output is `output.txt`.
Python 3 and its standard library suffice. All arithmetic is exact integer arithmetic.

`check.py` computes each coefficient in two independent ways: multiplication of
part-size factors, and the logarithmic-derivative recurrence for the complete
Euler product. Neither computation uses a conjectured congruence. It reproduces
all displayed OEIS S/T/U terms for each entry and checks 90 finite instances of
the separately stated Gauss congruences.

| OEIS entry | Source modulus and range | Prime cases through 200 | Result |
| --- | --- | ---: | --- |
| A270913 | p², all primes | 46 | All pass |
| A008485 | p², primes p ≥ 3 | 45 | All pass |
| A270919 | p², all primes | 46 | Only p = 2 fails: a(2) = 12 ≢ 6 (mod 4) |
| A008705 | p², primes p ≥ 3, both a(p) and a(2p) | 45 for each clause | All pass; maximum index 398 |
| A291697 | p³, primes p ≥ 5 | 44 | All pass |

The last prime in this range is 199. The excluded boundary cases A008485 at
p = 2 and A291697 at p = 3 are also checked and fail the corresponding formula.
`ALL_CHECKS_OK` asserts agreement with these explicitly stated outcomes; it
does not assert that every source conjecture is true, nor prove an infinite claim.

The definitions and conjectures are Peter Bala's OEIS comments of 2023-04-18:

- <https://oeis.org/A270913/internal>, revision 40, 2026-02-23.
- <https://oeis.org/A008485/internal>, revision 59, 2023-04-20.
- <https://oeis.org/A270919/internal>, revision 27, 2023-10-03.
- <https://oeis.org/A008705/internal>, revision 64, 2026-05-30.
- <https://oeis.org/A291697/internal>, revision 17, 2023-10-04.

The S/T/U terms embedded in the program are fidelity data from those entries.
OEIS terms and attribution information are at
<https://oeis.org/wiki/Legal_Documents>.
