# Exact checks for the minimal non-statistical manifold question

Source: Aliasghar Sarizadeh, *Minimal Oscillation of Cesàro Averages Implies Non-Statisticality*, arXiv:2605.31451v2, Section 1 (unnumbered question), and the displayed pushforward/integral identity in the proof of Theorem A, Section 2.

Run with Python 3 and SymPy:

```sh
python3 docs/reports/sarizadeh-2605.31451v2-minimal-lebesgue-historic/check.py
```

The output is stored in `output.txt`. The program uses integer and rational arithmetic, plus a symbolic rational-function normalization; it uses no floating-point values as evidence.

The fidelity check reproduces the source's finite-time pushforward/integral identity on every self-map of sets of one through four elements, with every atomic basis probability and one nonuniform full-support probability, at lengths 1 through 12. This identity does not require manifold assumptions. It is not a numerical test of dominated convergence or of the source's circle obstruction.

The proposed settlement route is checked through the exact empirical-measure conjugacy and endpoint-shift identities on all permutations of at most four states, at lengths 1 through 12; the normalization and shifted-weight bounds for the two-sided geometric mixture; and exact endpoint inequalities for a symbolic alternating-block model. Infinite geometric tails are included by their exact sums. The finite cyclic mixture check uses the sum of every congruence class in the infinite two-sided series, not a truncation.

These checks support algebraic identities in the proposed general proof. They do not establish a minimal non-uniquely ergodic torus homeomorphism, a nonatomic full-support measure concentrated on its historic set, measure transport by a homeomorphism, or the claimed manifold existence. Those remain mathematical proof obligations. The alternating-block model is a symbolic-shift illustration and is not a torus example.

The final output line is `ALL_CHECKS_OK`.
