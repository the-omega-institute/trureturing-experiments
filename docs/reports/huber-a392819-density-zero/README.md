# Huber A392819: exact finite checks

Source: Felix Huber, [OEIS A392819](https://oeis.org/A392819/internal), internal revision 23, March 2, 2026. The entry defines positive integers n for which p + v_p(n) divides n for each prime divisor p, and conjectures natural density zero. The 51 initial terms and two examples are attributed to that entry; OEIS data is available under [CC BY-SA 4.0](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement).

Run with Python 3 and no third-party dependencies:

```sh
python3 docs/reports/huber-a392819-density-zero/check.py
```

`check.out` is the program output. All arithmetic used for evidence is integer arithmetic or exact rational arithmetic with `fractions.Fraction`.

The check reproduces every source term through 2610 and the examples 54 and 840. Independent trial division checks the sieve factorization through 3000. Through 100000 it checks the defining condition against lcm divisibility and verifies that every prime divisor p of a term yields p^2 | n or p(p+1) | n. For cutoffs y = 2, 3, 5, 10, 20, 50 and every integer 1 <= X <= 100000, it checks the smooth-number bound and

c(X)/X <= 2^(number of primes below y) * floor(sqrt(X))/X + 2/(y-1).

At six endpoints per cutoff it also checks the exact divisor-union counting bound and rational reciprocal tail bound. The cases 12 and 54 exclude, respectively, replacing lcm divisibility by product divisibility and treating these numbers as a subset of refactorable numbers. Membership at 1 is vacuous; 0 is excluded.

These are finite consistency checks. They establish neither the infinite density limit nor a Lean proof. The proposed full argument uses a fixed smooth-number cutoff followed by a summable divisor tail.
