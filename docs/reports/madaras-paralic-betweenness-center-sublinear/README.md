# Tree betweenness centers: exact finite checks

Target: the final question of Problem 1 in Tomáš Madaras and Matúš Paralič,
*Betweenness centers of graphs*, arXiv:2609.09342v1: “Is f(n) = o(n)?”
The proposed general argument gives |C_b(T)| < 2√(2(n−1)) + 3 for n ≥ 2.
This program is supporting evidence, not a proof of that assertion.

Run with Python 3 and NetworkX 3.5:

```sh
python3 -m pip install networkx==3.5
python3 check.py > output.txt
```

`check.py` reproduces Example 1 and every row of Table 1 through order 13,
covering 2,288 nonisomorphic trees. It computes the literal unnormalized
betweenness scores by enumerating unordered vertex pairs and recovering their
unique paths, independently computes component sizes after vertex deletion,
and checks the cut identity at every vertex. It also checks the centroid
lower bound on the maximum score, path containment of all maximizing vertices,
the imbalance-spacing inequalities (including gaps in the center), and
reciprocal packing for every integer threshold from 1 to n−1.

All evidence uses integer arithmetic or `fractions.Fraction`. The radical bound
is tested by (k−3)² < 8(n−1) only when k > 3; the complementary cases satisfy
it directly. The singleton tree is checked separately. No floating-point
approximation is used. The final line of a successful run is `ALL_CHECKS_OK`.

The source definitions and table entries identify the fidelity target. The
algorithm and checks are independent implementations. Only the first 13 rows
are reproduced; no claim about exhaustive verification through order 20 or
about a universal theorem follows from this computation.
