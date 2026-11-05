# Unitary-orbit amplification: source fidelity and padding mechanism

Target: Marcoux–Sarkowicz–Zhang, arXiv:2508.13834v1, Question 4.16,
noncoprime clause (85). The expected answer is negative, using the
existence seed in Proposition 4.14 and scalar padding.

Run with Python 3 and SymPy:

    python3 check.py

All acceptance calculations use integers, fractions or exact SymPy algebra.
`output.txt` contains the output. `ALL_CHECKS_OK` means the listed checks
passed; it does not certify a counterexample to (85).

The program checks:

- The source's Section 3.1 ordered spectral matching and repetition formula,
  using an independent threshold matching algorithm (404 repetition cases).
- Amplification/padding permutation identities, including multiplicities
  2, 3 and 4, with padding sizes 1 through 4 (32 identities).
- Exact block residuals and polar-compression estimates on three rational
  unitary examples, with nonnormal inputs and resolvent bounds.
- The conditional family budget for every real q >= 1: seed gap at least
  1/q, lambda = 1+8q, padding loss 1/(8q^2), perturbation 1/(32q), and
  retained gap greater than 1/(4q). This is conditional on an actual seed.
- Scalar joint commutants for the irreducibility construction in dimensions
  2 through 8, and the gcd specialization on specified dimension windows.

## Unresolved seed check

No concrete seed dimension, entries or orbit distances are certified here.
The source's Proposition 4.14 and Marcoux–Zhang, JFA 280 (2021), 108778,
Theorem 5.1 are existence arguments. The latter explicitly says its proof
does not tell how to find the integers or matrices. The former leaves the
power of two unspecified and uses CAR algebra and dimension-drop embeddings.
The numerical request to compute delta_1 and delta_3 has therefore not been
fulfilled. Sampling a different pair would not check the published seed.

The report supports the conditional padding route, not a complete independent
check of the proposed counterexample family. A formal settlement needs a
proved seed, either through its source existence argument or through an
explicit exact witness; no new axiom or assumed seed is justified by this
report.

Sources: https://arxiv.org/abs/2508.13834v1 and
https://doi.org/10.1016/j.jfa.2020.108778 (publisher PDF available at
https://hdl.handle.net/10012/20320).
