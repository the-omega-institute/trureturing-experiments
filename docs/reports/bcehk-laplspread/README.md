# Laplacian spread of digraphs: exact checks

`check.py` checks definitions and candidate settlements of Conjectures 14 and 20
in Barrett–Cameron–Evans–Hall–Kempton, *On the Laplacian spread of digraphs*,
[arXiv:2206.15410v1](https://arxiv.org/abs/2206.15410v1), Sections 2–4 and 6.
It uses Python 3 and SymPy 1.14.0. Run `python3 check.py`; `output.txt` contains
the output. All certificates use integers, rational arithmetic and exact
polynomial root isolation. No floating point is used as evidence.

Fidelity checks reproduce the directed-cycle spreads in Corollary 3 for orders
2, 3, 4 and 6, and certify the spread of Figure 3 rounds to 5.035 using a
rational enclosure of the extremal root of its compressed characteristic
polynomial. The five vertices of Figure 3 are indexed left, upper middle,
lower middle, upper right, lower right.

For Conjecture 14, all 2, 10 and 152 labeled balanced digraphs of orders 2, 3
and 4 have spread zero or at least one. The degree Rayleigh identity is checked
for every such graph. The restricted-spectrum variance expression and its
strict lower bound are checked at 930 degree/order pairs with 5 ≤ n ≤ 64.
These finite checks do not establish the universal conjecture.

For Conjecture 20, the source's Figure 3 adjacency matrix has zero-sum integer
vectors x = (2,26,-26,-1,-1) and y = (0,0,0,1,-1), with Rayleigh values
3408/679 and 0. Replacing every vertex by an independent cluster of m vertices,
and every arc by all m² arcs between its clusters, gives an unweighted simple
digraph of order 5m. For m = 1,…,8, the script checks directly that the lifted
vectors have Rayleigh values m·3408/679 and 0. This certifies excess at least
13m/679 over the order in these cases. The preregistration's universal
blowup argument is ordinary mathematics and has not been checked in Lean.

The final output is `ALL_CHECKS_OK`. The checks are sanity and source-fidelity
evidence; neither candidate settlement is asserted to be kernel verified.
