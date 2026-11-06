# Binary parameters and star systems

Run `python3 check.py` with Python 3. The program uses only the standard library, integer graph searches and exact rational arithmetic. Expected output is in `check-output.txt`; the last line is `ALL_CHECKS_OK`.

The source is Contreras, Reyes and Villarreal, *Systems of parameters consisting of linear forms for monomial ideal quotients*, arXiv:2609.30081v2, Conjecture 4.18. The fidelity check reproduces Example 5.10's 11 vertices, 35 edges, independence number 2, chromatic number 6 and the full-rank restrictions of its displayed rational parameter matrix. Its printed expression for the covering number uses H where these counts concern G; that expression is not used as evidence here.

For B = M(C5), form the categorical graph product X = K3 × B, where (i,b) and (j,c) are adjacent exactly when i differs from j and bc is an edge of B. G is the complement of X. X has 33 vertices and 120 edges; G has 408 edges. Exhaustive DSATUR enumeration with the endpoints of an edge assigned colours 0 and 1 gives exactly one 3-colouring. Permuting colours gives precisely the six first-coordinate colourings.

Assign the column vectors (1,0), (0,1), (1,1) to the three layers. Every restriction to an independent pair of G has determinant 1 or -1. For every possible binary parameter colouring, each X-neighbourhood uses two column types; hence no row can have a centre coefficient 1 and vanish on that entire neighbourhood, as a star row requires. All ordered independent pairs of centres are thereby excluded. The same checks are run for independent-twin blow-ups with t = 1, 2, 3, 4. The local rigidity check examines all 27 fibre functions and their 66 compatible ordered pairs.

The rank-to-radical implication, dimension calculation and all-field infinite-family argument require the mathematical proofs stated in the preregistration. These finite computations are corroboration, not Lean kernel verification or a proof by testing every field or every t. The planned settlement is a symbolic product-rigidity argument for every connected triangle-free base graph with no 3-colouring. No third-party program source is used.
