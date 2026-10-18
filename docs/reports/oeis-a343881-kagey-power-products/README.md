# OEIS A343881: least m > k with k^x m^y an n-th power

Check for the open-problem preregistration refuting Kagey's conjecture in OEIS A343881: T(p, k) = A064549(k) for fixed k > 1 and sufficiently large p.

Run: `python3 check.py 100 43` (Python 3, standard library only; about 4 minutes).

The program decides "k^x m^y is an n-th power" exactly, from the prime exponents. It checks:
- the table rows n = 2, …, 11 for k = 1, …, 10 printed in the entry;
- T(p, 12) for the primes p < 48 (144 from p = 7 on), against A064549(12) = 72;
- that m = 72 works for no n ≤ 60 at k = 12;
- the corrected formula T(p, k) = k^(1 + 1/g), with g the gcd of the exponents of k, for 2 ≤ k ≤ 100 at the primes 31, 37, 41 and 43;
- the list of k ≤ 100 where the corrected formula differs from A064549(k), which starts at 12.

`run.log` holds the output of the run cited by the preregistration.
