# Same-family child cuts in the Han–Kazachkov model

The exact checker tests the conjecture in footnote 1 immediately before Lemma 4
of Han and Kazachkov, *The Strength of Root Cuts in an Extended Abstract
Branch-and-Cut Model*, arXiv:2610.31774v1.

Run with Python 3 and the standard library:

    python3 check.py

`check-output.txt` contains the output. The final line is `ALL_CHECKS_OK`.
`han_family_certificate.json` contains the affine Bellman action differences
and exhaustive-tree counts. No floating point arithmetic is used.

The checker reproduces the separation-placement example with
(c, ell, r, Z) = (1, 2, 4, 5), w(1) = 7/5 and a tail of 100. It checks the
source's two switching thresholds, including the cut-count tie-break.

The counterexample has cut strengths (61,81,95), branch gains (100,110),
and target 250. The unique optimal tree branches first, cuts with family 1
on both children, then cuts with family 3 on the left and family 2 on the right.
The two paths prove bounds 256 and 252. The rational base time is 6307/500.

An infinite family is parameterized by
41/20 <= a <= 43/20, 13/4 <= b <= 67/20,
1/1000 <= epsilon <= 1/100 and 1/1000 <= sigma <= 1/100.
The simple-time table, in history order, is

| History | Time |
|---|---|
| 0 | 1 |
| e1, e2, e3 | a, a+epsilon, a+2epsilon |
| 2e1, e1+e2 | a, a+epsilon |
| e1+e3, 2e2 | b, b |
| e2+e3, 2e3 | b+epsilon, b+2epsilon |
| total history n >= 3 | 100+n-3+(h2+2h3)*epsilon |

Every cut has separation time sigma. Finite boundary inequalities and the
symbolic tail increments establish the admissibility clauses. At the base
point epsilon=1/100, the tail agrees with the candidate's stated tail.

The Bellman checker evaluates 51 reachable states, selects a base policy,
and tests its affine action differences at all 16 parameter-box corners.
A separate generator enumerates all 770 trimmed proving trees with at most
12 nodes and at most two cuts on a path. It uses no Bellman policy or value
to prune these trees. Every competitor has strictly larger time at every
corner. Affineness extends this comparison to the whole box.

These bounds cover arbitrary finite trees: each node costs at least 1;
the displayed tree costs at most 257/20 < 13; any third cut has a child
costing at least 100; and continuing past a proved leaf strictly increases
time. Hence the displayed tree is unique and therefore cut-minimal.

The checker also verifies that all reachable prefix-threshold signs persist
under independent perturbations of the six gain/target parameters by at most
1/10. Lemma 4 is unaffected: after the common family-1 prefix, both children
continue with cuts from different families.

These are exact Python results and a mathematical outline, not a Lean proof.
The universal finite-certificate bridge and the tail/size reductions remain
obligations of the later Lean implementation.
