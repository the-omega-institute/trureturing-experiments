# Arbitrary-alphabet circuit lower bounds: exact finite checks

The target is Chen, Kjos-Hanssen, Koswara, Richter and Stephan, *Languages of Words of Low Automatic Complexity Are Hard to Compute*, arXiv:2510.07696v1, Question 5.1; conference DOI [10.4230/LIPIcs.FSTTCS.2025.24](https://doi.org/10.4230/LIPIcs.FSTTCS.2025.24), Question 47.

Run with Python 3, using only its standard library:

```sh
python3 docs/reports/chen-et-al-2510.07696-question-5.1-arbitrary-alphabet-circuits/check.py
```

The stored output is `output.txt`. A successful run exits 0 and ends with `ALL_CHECKS_OK`. All evidence uses integers or `fractions.Fraction`; there are no random samples or floating point calculations.

- **Source fidelity:** 411 instances of the source's Proposition 2.1, with nonempty base words of lengths 1–3, alphabet sizes 2–4 and thresholds 1/4, 1/3 and 2/5. The program checks the source's literal cyclic NFA, rejects any competing accepted word of the same length by divergence tracking, and checks the strict rational inequality. The source's empty-base-word endpoint is excluded from this check: the displayed cyclic construction has no state there and the empty word does not satisfy its strict threshold.
- **Grid support:** exhaustively enumerates all coefficient vectors in the printed ten degree/dimension/prime cases, including unreduced polynomials, all alphabet subsets of size at least two, zero variables, constants and polynomials vanishing on the whole grid. For a nonempty evaluation support it checks `support * r**degree >= r**n`. There are 3,812 polynomials and 25,612 polynomial/grid pairs, of which 25,051 have nonempty support.
- **Complement identity:** checks `1-P**(p-1)` is nonzero exactly where `P` is zero at every enumerated point. This is Chevalley's classical finite-field zero test.
- **Alphabet restriction:** computes exact NFA complexities by enumerating canonical accepting state paths and tests all binary words of lengths 0–4 in every ordered two-letter subalphabet of alphabets of sizes 2–4 (620 comparisons). An accepting path induces a transition subgraph and a singleton final set; every exactly accepting NFA has such a path. Thus this enumeration allows ambiguity but does not miss smaller NFAs.
- **Boolean encodings:** tests all 70 distinct ordered pairs of codewords of block widths 1–3; each bit restricted to the pair is a constant, the input bit or its negation. It also checks a specific NFA with two accepting paths, to distinguish exact from unique acceptance.

The candidate argument combines a general bounded-degree grid-support bound with the source's asymptotic density theorem (Corollary 2.10). These finite checks support the definitions and intermediate identities; they prove neither that density limit nor the unbounded circuit lower bound. The grid-support ingredient is classical: Alon–Füredi (1993), Theorem 5; see Bishnoi, Clark, Potukuchi and Schmitt, *On zeros of a polynomial in a finite grid*, Theorem 1.1, arXiv:1508.06020v2, DOI [10.1017/S0963548317000566](https://doi.org/10.1017/S0963548317000566).
