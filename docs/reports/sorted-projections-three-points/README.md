# Sorted projections of three-point multisets

The target is the fourth open-question bullet in Section 7 of Christopher
Gorham Lester and Raymond Hau Sørensen, “How many sorted projections determine
n points in k dimensions?”, arXiv:2611.06470v1: whether m_k(3)=2k for every k≥2.
Repeated points are allowed; one list of nonzero directions must work for all
triples simultaneously.

Run `python3 check.py` with Python 3 and SymPy. The stored output is
`output.txt`; a successful run ends in `ALL_CHECKS_OK`.

The program was independently written from the mathematical definitions in the
paper; the authors' ancillary program source was not inspected or reused.
Coordinates on the standard S3 representation are (a,b,-a-b). A matching tuple
produces the equations X d_j=P_j Y d_j. Within each parity class let N be its
size and n_P the multiplicity of P; put δ=max(0,max_P(2n_P-N)). On X=P Y for a
maximizer, at most 2k−δ independent conditions remain. Thus the full kernel has
dimension at least δ, and a rank lower bound 4k−δ forces containment in X=P Y.

The checks reproduce the paper's separating three-point lists in dimensions
2, 3 and 4, using all 216, 7776 and 279936 matching tuples after normalizing the
first permutation. Coordinate axes eliminate X. Rank is computed exactly over
the prime field with 65521 elements; hitting 2k−δ supplies a nonzero integer
minor and, with the preceding universal bound, the exact rational rank and
kernel containment. No floating-point calculations are used. If this modular
certificate fails, the program checks the rational kernel directly with SymPy;
no such fallback is needed for the reported source lists.

Additional checks cover a concrete instance of the lower-bound cancellation,
all 43065 multiplicity vectors of total 2k for 0≤k≤8, and rank witnesses for the
five residual cases. After pairing distinct same-parity permutations, each
parity leaves an empty class, a single permutation run, or one copy of each of
its three permutations. The residual cases are empty/run, run/run, odd run of
length ≥3/triad, run of length 1/triad, and triad/triad. The two exceptional
witness minors are 3 and −19440. The finite checks support a proposed uniform
argument; they do not prove the assertion for every ambient dimension.
