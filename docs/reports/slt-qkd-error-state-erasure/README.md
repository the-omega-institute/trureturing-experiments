# Public-flag error-state erasure in advantage distillation

`check.py` checks the operational conjecture in the paragraph immediately after
Theorem 2 of Stasiuk–Lütkenhaus–Tan, arXiv:2212.06975v3, Section III A.

Run with Python 3 and SymPy 1.14:

```sh
python3 docs/reports/slt-qkd-error-state-erasure/check.py
```

All evidence uses exact rational or symbolic algebra. The program checks the
three-qubit family's normalization and CHSH formula modulo its normalization
relations, then independently computes the conditional matrices and their spectra
at a rational nondegenerate witness. The public fair bit is retained in Eve's
four-dimensional input and output space. It also reproduces the source's Lemma 1
flag-swap identity for every message of length 1–3, and its comparisons
F² ≤ Q ≤ F, with F = 3/5 and Q = 9/25.

The witness uses a = 4/5, b = 3/5, s = 3/5 and rotation cosine/sine
9999/10001 and 200/10001. Each key outcome has positive probability and all four
conditional states are distinct. Trace-distance contraction and the triangle
inequality give a lower bound 38004/50005 on the sum of correct-state disturbance
and residual error-state distinguishability. Its positivity forbids simultaneous
exact preservation of both correct states and equality of the two error outputs.
The symbolic checks establish algebraic identities; the general channel
inequality is the standard mathematical implication documented above, not a
channel-search computation performed by this program.

`output.txt` is the program's complete standard output. The final line is
`ALL_CHECKS_OK`. These checks are not a Lean proof and do not establish a positive
key rate, refute Theorem 2, or exclude a replacement attack with the same observed
correlations. The general proof remains a subsequent formalization task.

Source: https://arxiv.org/abs/2212.06975v3 and
https://doi.org/10.1103/PhysRevA.111.022631.
