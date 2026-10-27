# OEIS A323206 parity: exact finite checks

Peter Luschny's [OEIS A323206](https://oeis.org/A323206), revision 52
(15 March 2026), defines the terminating hypergeometric array
A(n,k) = 2F1(-k,k+1;-k-1;n), for nonnegative integers n and k, and
conjectures that its entries are odd exactly when n is even or k+2 is
a positive power of two. The explicit conjecture retains the redundant
condition that n is odd in the second branch.

Run with Python 3 and its standard library:

```sh
python3 docs/reports/luschny-a323206-parity/check.py
```

The saved output is [check.out](check.out). The program uses integers and
`fractions.Fraction`; it performs no floating point computations.

The fidelity check reproduces all 72 displayed EXAMPLE entries with
0 <= n <= 8 and 0 <= k <= 7. Independently constructed Pochhammer products
agree with both published finite ballot formulas. The coefficient check
covers all 5151 coefficient cases 0 <= j <= k <= 100 and explicitly checks every
Pochhammer denominator for being nonzero. All 101 row sums equal the
shifted Catalan number, generated independently by first-return convolution
and checked against the central binomial formula. The conjectured parity
and the two residue classes are checked on all 2121 pairs with
0 <= n <= 20 and 0 <= k <= 100. The boundary cases n=0 and k=0 are included.
Successful execution exits 0 and ends with `ALL_CHECKS_OK`.

These are finite checks, not a proof for unbounded n and k. The fixed
example integers are transcribed from the cited OEIS revision; the source
is available under the [OEIS End-User License Agreement](https://oeis.org/LICENSE).
The check program is an independent implementation of the displayed formulas.
