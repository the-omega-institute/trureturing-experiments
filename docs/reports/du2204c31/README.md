# Du–You–Liu–Huang, maximum Sombor coindex of two-trees

Run `python3 check.py` with Python 3.12 or later. The committed output is in `check.out`.

The check generates two-trees by attachment on existing edges and reproduces the coindices of the two graphs in Proposition 2.1 at orders 5 through 16. It verifies the construction and directly enumerates unordered nonedges for both graphs at five orders, starting at 1,024. The comparisons use integer lower and upper enclosures of each square root. The check also verifies the arithmetic inequalities used in the uniform comparison for every integer parameter at least 32.

The target is the structural coindex maximizer from the conclusion of [arXiv:2204.02746v1](https://arxiv.org/abs/2204.02746) and [Conjecture 3.1 of the journal article](https://doi.org/10.3934/math.2023967).
