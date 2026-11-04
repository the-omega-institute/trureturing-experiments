# Cayley-code direct products

The check concerns Arumugam–Praeger–Rademacher, *Structure of Cayley Codes*,
arXiv:2606.28611v1, Proposition 6.2, Lemma 7.4, Proposition 7.5(b), and Question 4.

Run `python3 check.py` with Python 3.8 or later. The program uses only the standard
library and exact arithmetic modulo the primes 2, 3, and 5. `output.txt` is its
complete output; success requires exit 0 and final line `ALL_CHECKS_OK`.

The fidelity screen reproduces 114 cycle-code dimensions from Proposition 6.2 and
the source's zero/full local-code cases. The product screen computes codes directly
on unordered edges from local parity-check equations, verifies the two-element
fibers of the source edge map, and checks every basis generator of the subcode
constructed in Proposition 7.5(b) against the product constraints.

Independently, arc reversal and the simultaneous local/reversed local constraints
compute `r_G = rank((R_G - I)|D_G)`, with `D_G = U_G ∩ R_G(U_G)`. The program checks
`dim C_delta = dim C_G * dim C_H + r_G * r_H`, and equality of the constructed
subcode with the whole product code exactly when one reversal rank vanishes.

There are 606 product comparisons: 500 exhaustive choices of degree-two local
subspaces on ordered C3/C4 pairs; 90 deterministic further cycle pairs; six
complete-graph augmentation cases; four degree-four circulant augmentation cases;
and six nondegenerate equality examples pairing odd repetition cycles with a
degree-four augmentation circulant. The latter two classes support the proposed
infinite families. Exhaustiveness is asserted only for the stated C3/C4 scope.

These finite checks are supporting evidence for a general linear-algebra argument.
They establish neither a universal identity nor a Lean proof. Extension fields,
nonabelian groups, and unbounded orders are outside the computational screen.
