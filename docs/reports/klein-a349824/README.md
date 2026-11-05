# Klein's A349824 iteration

For a natural number n, F(n) is the number of prime factors of n times their sum, counting multiplicity in both quantities. F(0) = F(1) = 0. Mike Klein's conjecture says that every orbit reaches 0, a prime, 27, 30, or the cycle 28 ↔ 33.

Source: [OEIS A349824](https://oeis.org/A349824), revision 44. The source supplies 74 DATA terms and example trajectories from 4 and 6. Source attribution: N. J. A. Sloane; conjecture: Mike Klein; see the [OEIS end-user licence](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement).

Run with Python 3, without external packages:

```sh
python3 docs/reports/klein-a349824/check.py
```

The program uses integer arithmetic throughout. A sieve computes the multiplicity count and sum; independent trial division checks its values for 0 through 10,000. All 74 source terms and both trajectories agree. The complete F(0), …, F(72) table stays in that range and reaches the source target set within 11 steps.

The program checks one-step descent for every composite n from 73 through 200,000 except twice-primes; F(2q) ≤ q + 15 for every odd composite q from 39 through 200,000; and descent within three steps for 2p for every prime p from 37 through 200,000. Every orbit starting from 0 through 200,000 reaches the target set; the maximum observed number of steps to the set is 22, first at 172803. The output is in `check.out` and ends with `ALL_CHECKS_OK`.

These are finite exact checks. They do not prove convergence for all natural numbers. The proposed universal argument uses product versus sum estimates on arbitrary prime-factor lists, the modulo-3 obstruction to three consecutive odd primes at distance two, and strong induction. The finite table supplies only the base case.
