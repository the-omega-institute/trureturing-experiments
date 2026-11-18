# Exact checks for Hasni et al., Conjecture 1

Source: R. Hasni, F. Movahedi, H. Kamarulhaili and M. H. Akhbari,
“Sharp Upper Bound for Harmonic Index of Trees with Given Total Domination Number”,
Ars Combinatoria 159 (2024), 179–186, Conjecture 1, equation (3).
DOI: <https://doi.org/10.61091/ars159-15>.

Run with Python 3.12 or newer and SymPy 1.14:

```sh
python3 docs/reports/hasni24c1/check.py
```

`check.out` is its complete standard output. All arithmetic uses integers,
`fractions.Fraction` or exact SymPy expressions.

The fidelity checks reproduce source Lemma 1 on P4, P8, P12 and P16,
and the stated P4/S4 base-case harmonic indices. The program constructs
adjacency independently, tests connectedness and the tree edge count,
and exhaustively minimizes total dominating sets. It checks the proposed
family B_r for 1 ≤ r ≤ 32 and the source's equality template for
1 ≤ r ≤ 4 and 1 ≤ s ≤ 3.

For B_r, a hub with r leaves and one arm of length four, it checks
γ_t = 4 and the exact harmonic sum. Exact symbolic algebra verifies

\[
F(r+5,4)-H(B_r)=\frac1{12}-\frac1{(r+2)(r+3)},
\]

and

\[
F(r+5,4)-H(B_r)-\frac1{20}
=\frac{(r-3)(r+8)}{30(r+2)(r+3)}.
\]

After substituting r = k + 3, numerator and denominator have nonnegative
coefficients and the denominator has positive constant coefficient.
The output distinguishes the finite graph checks from these symbolic
arithmetic certificates. No Lean settlement proof or literature novelty
claim is part of this experiment.
