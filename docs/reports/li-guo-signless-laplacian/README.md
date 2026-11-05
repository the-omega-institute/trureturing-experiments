# Li–Guo signless Laplacian inequality

Source: Wen-Jun Li and Ji-Ming Guo, “A new conjecture on Laplacian eigenvalues and degree sequence of a graph”, Filomat 39:13 (2025), 4551–4564, Conjecture 2.22. DOI: https://doi.org/10.2298/FIL2513551L.

Run with Python 3 and SymPy 1.14.0:

    python3 docs/reports/li-guo-signless-laplacian/check.py

The saved output is `check.out`. Arithmetic uses integers, rational numbers and exact reduction modulo cyclotomic polynomials. Real-root isolation gives rational enclosures for the square sums of the positive and negative eigenvalues of the actual integer adjacency matrix. No floating-point value decides an inequality.

The check reproduces Lemma 2.5(1) for complete graphs of orders 2–8 and checks the cycle spectrum in Lemma 2.5(4) against direct characteristic polynomials for orders 3–35. It also reproduces the ordinary-Laplacian cycle inequalities in Theorem 2.14 and Remark 2 on that range.

For every cycle order 3–35 it checks the proposed signless truncated-square formula, including strict violations at all nine orders congruent to 3 modulo 4. In particular, the triangle has signless spectrum 4, 1, 1, so its truncated-square sum is 4 and the proposed upper bound is 3. These are finite exact checks; an all-orders assertion and the classification of regular graphs require mathematical proofs. No Lean proof is included.
