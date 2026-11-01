# A378145 floor-weighted Catalan row sum

Werner Schulte, [OEIS A378145](https://oeis.org/A378145), revision 11
(author attribution 2024-11-17; revision 2024-12-08), second FORMULA conjecture:

> Sum_{k=0..n} (1 + floor(k/2)) * T(n, k) = A000108(n+1).

The triangle is the Riordan array (1 + x C(x), x C(x)), where C is the
Catalan generating function. Its entries are coefficients of
(1 + x C(x)) (x C(x))^k. The source's separate binomial-shift conjecture
is outside this report's target.

Run with Python 3, without external packages:

```sh
python3 docs/reports/a378145-floor-weighted-catalan/check.py
```

The program uses only integers and exact fractions. It independently builds
Catalan numbers by convolution, checks their binomial formula, computes
triangle entries from the source's piecewise binomial formula, and compares
them with truncated formal-series multiplication.

The stored output in check.out records:

- all 66 published DATA terms reproduced;
- 3,321 entries through row 80 matching the Riordan coefficient definition;
- 3,240 above-diagonal coefficients vanishing;
- 3,160 instances of the source recurrence verified;
- all 81 floor-weighted row sums for 0 <= n <= 80 equal to Catalan(n+1).

These are finite checks. They do not establish the universal conjecture.
The diagonal formula branch supplies T(0,0)=1; every denominator used in
the off-diagonal branch is checked positive. The final output line is
ALL_CHECKS_OK, and the program exits 0.

Source data attribution: Werner Schulte and the OEIS Foundation Inc.;
[OEIS End-User License Agreement](https://oeis.org/LICENSE).
The check program is an independently written implementation.
