# Alternating spanning-tree recurrences

Exact finite controls for Paul Raff, *Spanning Trees in Grid Graphs*, arXiv:0809.2551v1, Conjecture 2. The expected settlement is a general proof; these checks do not prove the universal statement or establish literature novelty.

Run with Python 3 and SymPy 1.14.0:

```sh
python3 docs/reports/grid-spanning-tree-alternating-recurrence/check.py
```

The output is stored in `output.txt`. A successful run exits 0 and ends with `ALL_CHECKS_OK`.

The graph product is Cartesian. A separate integer Bareiss determinant computes a Laplacian cofactor of the actual product graph. Every Bareiss division is checked for exactness. The spectral formula uses exact algebraic Laplacian eigenvalues and the polynomial recurrence S_0=0, S_1=1, S_(n+1)=(lambda+2)S_n-S_(n-1). Resultants construct the full signed-product annihilator without choosing square roots. Integer coefficients, strict alternating signs and positive roots are checked exactly, with Sturm counting for the last property.

- **Fidelity:** Raff Theorem 6.3 gives the ladder recurrence T_n=4T_(n-1)-T_(n-2), with T_1=1 and T_2=4. Independent cofactors reproduce 1, 4, 15, 56, 209 and check it through index 19. The P_3 and P_4 initial terms and the P_4 order-eight coefficients also reproduce the source's table.
- **All small graphs:** all 76 labelled simple graphs on 0 through 4 vertices are included (44 connected; 32 empty or disconnected). For each, six product-graph cofactors, at path lengths 1 through 6, agree with the exact spectral formula. The formula-generated values through index 19 satisfy the alternating integer annihilator; this longer recurrence check is not an independent cofactor check for every graph.
- **Extra controls:** for P_4 and K_(1,3), independent cofactors additionally check lengths 7 through 19. P_4 has irrational Laplacian eigenvalues; K_(1,3) has repeated eigenvalues and coincident modes. Forty-two further cofactors check K_k times P_n, with 2<=k<=8 and 1<=n<=6, against k^(k-2) S_n(k)^(k-1).

All evidence uses Python integers and exact SymPy expressions; no floating point is used. There are 543 cofactor evaluations, 524 cofactor/formula comparisons and 1089 recurrence equations in this run. The proposed general argument still requires the matrix-tree/counting bridge, path determinant, algebraic-integrality descent and exponential-mode annihilation. No Lean proof is supplied here.

Sources: <https://arxiv.org/abs/0809.2551v1>; the repeated conjecture in Raff's dissertation, <https://doi.org/10.7282/T3T72HM0>. The script is original checking code; the papers are referenced, not redistributed.
