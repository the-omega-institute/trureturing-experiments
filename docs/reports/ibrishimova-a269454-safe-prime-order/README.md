# Ibrishimova's safe-prime product-order conjecture

Source: [OEIS A269454](https://oeis.org/A269454/internal), revision 54, second comment. A safe prime is a prime p for which (p−1)/2 is prime ([A005385](https://oeis.org/A005385)); hence 5 is included and 3 is excluded. The conjecture says that distinct safe primes p and q, neither congruent to 7 modulo 8, satisfy ord_(pq)(2) = phi(pq)/2.

Run with Python 3 and its standard library:

```sh
python3 docs/reports/ibrishimova-a269454-safe-prime-order/check.py
```

The program uses integer arithmetic exclusively. It reproduces every displayed A269454 term through 7523 and checks the source's universal half-totient exponent and upper bound for all 300 distinct odd-prime pairs at most 101. Totients in that check are also counted directly using gcd.

For every unordered distinct pair from the 53 admitted safe primes at most 10000, it computes the order by reducing Euler's exponent and checks the conjectured equality. Each returned exponent is certified to give residue 1 and to have no quotient by a prime divisor giving residue 1. This checks 1378 pairs, including 52 containing 5. All 15 admitted pairs at most 200 are also checked by independent successive-residue enumeration. The excluded pair (7,23) illustrates why an unrestricted equality would fail.

`check.out` is the complete output. These are finite checks supporting a preregistration; they do not prove the universal conjecture. There is no Lean proof in this report.
