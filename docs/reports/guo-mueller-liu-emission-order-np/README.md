# Gaussian emission order and graph pathwidth

Run `python3 check.py` with Python 3. The program uses only the standard library,
integer arithmetic and `fractions.Fraction`. The stored output ends in
`ALL_CHECKS_OK`; assertions stop execution on any failed equality.

The source is Fucheng Guo, Frank Mueller and Yuan Liu, *Memory-Optimal Sequential
Synthesis of Multimode Gaussian Transformations*, arXiv:2609.05250v1,
https://arxiv.org/abs/2609.05250.

The fidelity check reconstructs the source's nearest-neighbour two-mode-squeezing
chain using exact rational values cosh(r)=5/3 and sinh(r)=4/3. It reproduces the
five-mode forward and reverse occupation profiles (2,2,2,2,1) and (5,4,3,2,1).
It also verifies Corollary 1's rank/count equality at every cut of every order for
all graph shears with at most four modes, and for additional non-shear symplectic
matrices. Ranks use exact Gaussian elimination over rational numbers.

The candidate check covers all 1,099 labelled simple graphs on one through five
vertices and all 124,469 vertex orders. It verifies symplecticity, the exact
closed-neighbourhood block supports, and the per-order identity
C(S_G,pi)=1+VS(G,reverse(pi)). The order-derived bags are checked for vertex and
edge coverage and the interval property. A separate exhaustive search over
introduce/forget bag sequences computes pathwidth directly and checks each
minimum against 1+pathwidth(G).

These finite checks support the proposed universal graph/shear correspondence.
They do not prove NP-hardness, certify a polynomial-time machine reduction, or
settle the unbounded claim. The preregistration's main route is a general
argument; enumeration is supplementary evidence.
