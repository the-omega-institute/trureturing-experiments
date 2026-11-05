# Krasniqi boundary-contact uniqueness: exact preregistration checks

Run with Python 3 and its standard library:

```sh
python3 check.py
```

The program uses only integer and rational arithmetic. Its expected output is
stored in `output.txt`, ending with `ALL_CHECKS_OK`.

The fidelity check reconstructs the coefficient polynomials in Krasniqi,
*On the Alzer–Berg problem: an optimal Bernstein boundary and a uniqueness
conjecture*, arXiv:2609.28707v1, equations (3.17)–(3.19). It independently
compares the recurrence through q9 with the truncated product of exponential
series, and checks the four displayed coefficient formulas symbolically.
The binomial transform constructs d0 through d8, including the prescribed
constant and linear density terms.

The expected-settlement check independently implements sparse polynomial
arithmetic, affine substitution and conversion to tensor Bernstein
coefficients. Its partition is the 238-leaf Appendix A of
`docs/develop/theory/KrasniqiBoundaryContact.md` at main-repository revision
458aaeee2237c03eff58e0bb81b72d1d7a7431a1. It reconstructs P6 from the source
recurrence; checks a full closed-face binary partition; verifies all 134
positive leaves and 104 exclusions; bounds the logarithmic relaxation with
an exact positive series and geometric tail; and checks the derivative at
t=16 on the whole parameter rectangle. The resulting rational lower
margins are 6583/11943936 and 45599851/571536000.

This checks exact algebraic premises for the universal analytic proof route,
not a finite sample of conjectural contact points. It does not prove
Conjecture 23 or identify the exact boundary numerically. The source-faithful
improper integral, differentiated moment identity, boundary attainment and
necessary inequalities, positive filtered derivatives, Taylor remainder,
and monotonicity argument still need their analytic and Lean bridges.
The certificate is a helper for that general argument, rather than the
positive enumeration of a finite family as the main theorem.
