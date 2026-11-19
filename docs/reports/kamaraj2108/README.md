# Exact checks for Kamaraj, Johnson and Satheesh

Run `python3 check.py` with Python 3.12+ and SymPy 1.14.0. The reference output is in [check.out](check.out).

Example 3.2 checks the source conventions by reproducing the reported weighted inverses and the failure of reverse order there. The symbolic rank-one family checks the reverse-order identity, both range failures and the two block ranks for each nonzero real parameter. The final check verifies the Gram determinant identity used in the dimension-two argument. [arXiv:2108.05873v1](https://arxiv.org/abs/2108.05873v1) · [Filomat 37(3), 699–709](https://doi.org/10.2298/FIL2303699K)
