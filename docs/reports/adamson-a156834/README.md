# Adamson's A156834 prime-equality conjecture

For positive n, OEIS A156834 defines

a(n) = sum over positive divisors d of n of phi(d) binomial(d+n/d-2,d-1).

Its comment conjectures that, for n>1, a(n)=n if and only if n is prime.
The expected settlement is a general proof by nonnegative divisor contributions.
This report contains finite checks; it does not certify the unbounded conjecture.

Run with Python 3 and its standard library:

```sh
python3 docs/reports/adamson-a156834/check.py
```

The deterministic output is stored in `check.out`. A successful run exits 0
and ends in `ALL_CHECKS_OK`. All calculations use exact Python integers.

The checks reproduce all 51 terms printed in OEIS revision 20 and the first
11 triangle rows in the paper's section 5.2. The triangle is independently
constructed by Pascal addition and column interleaving through row 256;
every entry agrees with Proposition 6, and matrix multiplication by totients
agrees with the divisor-sum definition. Totients are checked against trial
factorization through 10000 and direct coprime counting through 200.

For every n from 1 through 10000, the program checks the classical totient
sum, the binomial equality criterion, and the integer excess identity.
For n from 2 through 10000, equality occurs exactly at the 1229 primes;
the 8770 composites give strictly positive excess. The excluded boundary
n=1 has a(1)=1 although 1 is not prime.

Sources and attribution:

- Gary W. Adamson, [OEIS A156834](https://oeis.org/A156834/internal),
  revision 20, submitted February 16, 2009. The explicit divisor sum is
  credited to Seiichi Manyama, April 22, 2021. The numerical prefix is
  transcribed from the entry. See the [OEIS license agreements](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement).
- el Houcein el Abdalaoui, Mohamed Dahmoune, Djelloul Ziadi,
  [On the transition reduction problem for finite automata](https://arxiv.org/abs/1301.3751v1),
  arXiv:1301.3751v1, January 16, 2013, section 5.2, Proposition 6.
  The triangle's eleven displayed rows supply the second fidelity target.

The scripts and computed output contain no Lean proof or refutation.
