# Bala's Tribonacci prime-power congruence

Run `python3 check.py` using Python 3. The program uses only standard-library
integer arithmetic and exits 0 with final line `ALL_CHECKS_OK`. Its captured
output is `output.txt`.

The source is Peter Bala's unnumbered FORMULA conjecture in
[OEIS A000073](https://oeis.org/A000073), Dec 28 2022, current revision #742;
[OEIS A106282](https://oeis.org/A106282), revision #18, defines the prime set.
The sequence starts `0,0,1`, not `0,1,1`.

Source fidelity checks reproduce all 38 displayed A000073 terms and all 37
A106282 primes below 620. The conjecture grid covers every eligible prime
below 620, every `n=1,...,20`, and every `k=1,...,6` (4440 cases). An additional
888 cases include `n=p,p^2,p^3` and compare independent polynomial-quotient and
companion-matrix exponentiation. These two algorithms also agree with direct
recurrence on 2107 small index/modulus pairs. A separate direct-recurrence
check covers 30 conjecture cases for `p=3,5`. Negative controls at ineligible
primes 2 and 7 confirm that the prime restriction matters.

All checks use exact integers. They provide finite evidence and source
fidelity, without establishing the universally quantified conjecture.
