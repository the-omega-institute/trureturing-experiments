# Mathar's seven-term recurrence for OEIS A124753

`check.py` checks the rational sequence

\[
a(3m+s)=\frac{(s+1)\binom{4m+s}{m}}{3m+s+1},\qquad m\ge0,\quad 0\le s<3,
\]

against R. J. Mathar's seven-term recurrence conjecture in the FORMULA field of
[OEIS A124753, revision #71](https://oeis.org/search?q=id:A124753&fmt=text).
The formula is attributed to Paul Barry and the recurrence to R. J. Mathar.
OEIS data and formulas are attributed under the [OEIS End-User License Agreement](https://oeis.org/LICENSE).

Run with Python 3 and SymPy:

```sh
python3 docs/reports/mathar-a124753-rec/check.py
```

`output.txt` is the program's exact output. The program verifies all 39 published
data terms, the entry's convolution formula for indices 1 through 150, and the
conjectured recurrence for indices 6 through 300 using `fractions.Fraction`.

For each residue `s=0,1,2`, write `n=3m+s`, and for each shift `j=0,...,6` put
`t=(s-j)//3`, `u=s-j-3t`. The common factorial scale is

\[
C_s(m)=\frac{(4m+s-8)!}{m!(3m+s+1)!}.
\]

With falling products \(F_r(z)=\prod_{h=0}^{r-1}(z-h)\) and rising products
\(R_r(z)=\prod_{h=1}^{r}(z+h)\), the multiplier is

\[
P_{s,j}(m)=(u+1)F_{-t}(m)R_{8-j+t}(4m+s-8)F_j(3m+s+1).
\]

The program checks the factorial bridge `a(3m+s-j)=C_s(m)*P_{s,j}(m)` for all
21 shifts and `m=2,...,50`. It verifies the three degree-at-most-11 polynomial
identities obtained by weighting these multipliers by the conjectured
coefficient polynomials. Both exact SymPy expansion and a separate polynomial
implementation using only Python integer addition and convolution give twelve
zero coefficients for each identity.

The polynomial certificates are exact identities. The sampled sequence and
factorial-bridge checks are finite diagnostics. They are not a Lean proof of
the universally quantified recurrence; that proof is a separate deliverable.
