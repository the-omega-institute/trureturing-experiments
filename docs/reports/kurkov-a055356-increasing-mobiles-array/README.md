# Kurkov's A055356 increasing-mobiles array

Run `python3 check.py` with Python 3; only the standard library is used. All arithmetic is exact integer arithmetic. The stored output is `check.out`, whose final line is `ALL_CHECKS_OK`.

The check reproduces all seven displayed rows of OEIS A055356 (including terminal zeros), independently counts increasing mobiles with 2 through 8 vertices by increasing parent functions weighted by cyclic child orders, and reproduces the third derivative of Ma's mobile grammar. This derivative corresponds to OEIS row 4; Ma's printed row index is shifted by one relative to OEIS.

For 0 <= n <= 10 and 0 <= k <= 8, three independently implemented constructions agree: Kurkov's nonlinear binomial-convolution recurrence, direct iteration of the derivation D(a)=ab+xa² and D(b)=xab on integer multivariate polynomials, and the proposed linear differential recurrence. The monomial support invariant is checked before evaluation. The first column agrees with OEIS rows n+2 for every checked n.

These are finite checks, not a proof of the universal identity. No Lean proof or admission is claimed.

Sources:

- [OEIS A055356](https://oeis.org/A055356), revision 45: row data, Whitmore's row recurrence, Fischer's coefficient recurrence, and Kurkov's final FORMULA conjecture. OEIS content is subject to the [OEIS end-user license](https://oeis.org/LICENSE).
- Shi-Mei Ma, [Some combinatorial sequences associated with context-free grammars](https://arxiv.org/abs/1208.3104v2), Theorem 6(c12): the grammar used for the indexing check.

The check program is original code; no source-paper text or third-party code is included.
