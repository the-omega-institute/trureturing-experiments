# Alternating repunit multipliers: OEIS A391473

`check.py` performs exact finite checks of Gonzalo Martínez's two conjectured formulas in [OEIS A391473](https://oeis.org/A391473), revision 28, and the alternating-integer definition in [A110305](https://oeis.org/A110305), revision 15. A391473 attributes the sequence to Martínez, Dec 10 2025; its current revision timestamp is 2026-07-02T01:02:11-04:00. A110305 attributes its definition to Walter Nissen. The quoted sequence data and definition are OEIS material under [CC BY-SA 4.0](https://oeis.org/LICENSE).

Run with Python 3.12 or newer, using no external packages:

```sh
python3 docs/reports/martinez-a391473-alternating-repunit-multipliers/check.py
```

The stored output is `output.txt`. The program exits 0 and ends with `ALL_CHECKS_OK` when all assertions pass.

- Increasing-multiplier brute force reproduces the first eight published terms: 1, 11, 19, 371, 1891, 36911, 189091, 3690911. It does not use the conjectured formulas as a stopping bound.
- A separate breadth-first search on remainder and last-digit parity reproduces the least alternating multiples for repunit lengths 1 through 6. The first prefix reaching a state has the least length and then the least lexicographic value; every continuation depends only on that state.
- Both formulas agree with all published terms in their domain, lengths 3 through 17.
- For every length 3 through 200, the explicit high and low blocks have the claimed digits, their concatenation equals the proposed multiplier times the repunit, the product alternates and has exactly twice the length minus two digits, and the geometric-sum identities hold. For even lengths the paired-digit parities and individual carry equations are checked as well.

All evidence uses Python integers and `fractions.Fraction`; no floating point, numerical tolerance, network query, random sampling, or third-party program source is used. These finite results do not prove universal minimality. The proposed settlement requires an argument for arbitrary length showing that every smaller positive multiplier fails.
