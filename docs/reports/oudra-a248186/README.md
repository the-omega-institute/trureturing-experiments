# Oudra's cubic floor conjecture (OEIS A248186)

Source: [OEIS A248186, revision 18](https://oeis.org/A248186/internal).
Clark Kimberling defines a(n) as the least positive integer k such that

\[
\frac1{18}-\sum_{h=1}^{k}\frac1{h(h+1)(h+2)(h+3)}<\frac1{n^3}.
\]

Ridouane Oudra's FORMULA conjecture asserts, for every integer n > 2,

\[
a(n)=\left\lfloor rn+\frac{r^2}{n}-1\right\rfloor,
\qquad r=3^{-1/3}.
\]

Run from the experiments repository root:

```sh
python3 docs/reports/oudra-a248186/check.py
```

The program requires only Python 3's standard library. Its integer and
`fractions.Fraction` calculations use no floating point. It constructs the
finite sums directly, checks the telescoping residual through k = 10,000,
reproduces the source's a(5) = 2 example, its five maximum-index table rows,
and all 68 displayed DATA terms. It checks the conjecture for all
3 ≤ n ≤ 10,000. Integer cube comparisons isolate r between rational
endpoints of width 1/10^80; the formula has the same floor at both endpoints
for every tested n. Thus its floor values are certified exactly.

`output.txt` contains the deterministic output, ending in `ALL_CHECKS_OK`.
These are finite checks supporting a proposed general argument; they do
not establish the universal assertion or supply a Lean proof.

The source's NAME omits the word “positive” for k. Its DATA starts with
a(1) = a(2) = 1, and its Mathematica program searches `Range[z]`, so the
check searches positive k. The conjecture is tested only on its stated
domain n > 2.
