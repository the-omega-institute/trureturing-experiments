# Exact-rate independent-set check

This report checks the construction nominated for Koch-Hyde–Olive, arXiv:2609.00382v2, Section 5, Problem 5.1. The Python program uses signed integer letters and a stack implementation of free reduction, so every equality test is exact.

It checks prefix-freeness of the conjugator family, the source lower count and exact conjugate lengths for ranks 2 and 3, the common-prefix junction endpoints, bounded nontrivial products of the resulting generators, and the elementary all-words upper bound with `fractions.Fraction`. The target rates are recorded with SymPy’s exact `sqrt`.

Run:

```text
python3 check.py
```

The captured output is `check.out`; the structured part is `result.json`. The final line is `ALL_CHECKS_OK`.
