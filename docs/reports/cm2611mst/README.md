# Exact checks for Chambrion–Mendico Conjectures 1 and 2

Source: Thomas Chambrion and Cristian Mendico, *Numerical investigation of a MST conjecture for random Haar matrices*, arXiv:2611.14091v1, Sections 1.2–1.3 and 2.1–2.2.

Run `python3 docs/reports/cm2611mst/check.py` from the experiments repository root. The program uses only the Python standard library and exact integers, rational numbers and finite rational linear combinations of square roots. Its output is in `output.txt`.

The fidelity checks reproduce the source's exact off-diagonal second moment, d v_d/(d²−1), at every dimension 2 through 12, and its limiting spectral-variance coefficient 9/100. Actual rational orthogonal conjugations of the exact spectrum k^(3/2) verify the column-radius identity, for every column at dimensions 2 through 12. The leading constant is checked by the exact equality (3/10)² = 9/100 and its reciprocal 10/3.

The proof-route checks verify the grouped Dirichlet inverse-moment bound at every dimension 16 through 1024. They enumerate all 1783 complete weighted graphs with 2 through 4 vertices and edge weights in {1,2,3}, or 5 vertices and weights in {1,2}. For every graph and every nonempty proper anchor set {0,…,k−1}, a fixed anchor chain plus the cheapest attachment of each remaining vertex is a tree bounding the exact minimum spanning-tree weight from above. Counting edges cheaper than each of four rational thresholds bounds that minimum from below. The exact minimum is independently obtained by enumerating all Prüfer trees.

These finite checks support the algebra and the two deterministic tree inequalities. They prove neither the asymptotic Haar statement nor the conditional distribution and concentration lemmas. There is no floating-point evidence or random-matrix simulation. The positive settlement is planned through a general analytic argument, with these computations serving only as diagnostic checks.
