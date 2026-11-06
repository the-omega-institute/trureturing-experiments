# Minimal transitive factorizations on K_{1,1,n-2}

Source: Cordelia Yuqiao Li and Ricky Ini Liu, *Minimal transitive factorizations supported on quasi-threshold graphs*, [arXiv:2609.31982v1](https://arxiv.org/abs/2609.31982v1), §5.4, second unnumbered conjectural formula.

Run `python3 check.py`. The program requires only Python 3's standard library and uses integers and `fractions.Fraction` throughout. Its deterministic output is in `output.txt`; success ends with `ALL_CHECKS_OK`.

The permutation-word dynamic program records both the product and the exact set of used edges. It counts identity words of length 2n−2 whose spanning support is connected. Support may be a proper subgraph of the allowed graph. As source fidelity checks, it reproduces the §4.2 example for K4 (16 factorization trees, total weight 96, count 2880) and the §2.2 star counts for n=2,…,6.

For K_{1,1,n−2}, independent word counts at n=3,…,6 are 24, 900, 47040 and 3175200. A separate exhaustive weighted-tree enumeration agrees with the proposed classification and factorial formula. Every edge is oriented from its smaller label to its larger label; prefix flow inequalities give the exhaustive bound 1≤w≤n−1 for every weight. The n=5 count requiring every allowed edge to occur is 1920, also agreeing with Shapiro, arXiv:2610.14373v1, Theorem 6; that count differs from the conjecture's 47040.

The scalar double sum, Vandermonde reduction, factorial value and individual telescoping identities are checked for m=0,…,100. These finite checks support a preregistration draft. They establish neither the unbounded claim nor a formal bridge from the original factorization count to weighted trees.
