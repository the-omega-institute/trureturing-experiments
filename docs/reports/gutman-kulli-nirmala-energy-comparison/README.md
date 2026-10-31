# Exact Nirmala energy comparison checks

Gutman and Kulli, “Nirmala energy”, Open Journal of Discrete Applied Mathematics
4(2) (2021), 11–16, DOI: [10.30538/psrp-odam2021.0055](https://doi.org/10.30538/psrp-odam2021.0055).
The matrices use the source's equation (5): an edge ij has weight
`sqrt(degree(i) + degree(j))`; all other entries are zero. Energy is the sum of
absolute eigenvalues with algebraic multiplicity.

Run with Python 3 and SymPy 1.14.0:

```sh
python3 docs/reports/gutman-kulli-nirmala-energy-comparison/check.py
```

The program enumerates all 76 labelled simple graphs on zero through four
vertices and checks 19 representatives under all vertex permutations. It verifies
Lemma 1's trace identities, Theorem 4 on every tree in that range, the exact
regular-graph scaling identity, and the expected strict energy comparison for
every graph with an edge. It also checks the stronger lower bound
`E_N(G) >= sqrt(2 delta_positive) E(G)`, where `delta_positive` is the least positive
degree. Proportional matrices certify sharp equality; in all other cases exact
rational root isolation certifies a strict inequality.

Only integer, rational and symbolic algebraic arithmetic enters the acceptance
criteria. The characteristic polynomial becomes rational after a verified exact
rescaling. SymPy's Sturm root intervals then enclose the energy; every root and
its multiplicity is counted. Edgeless graphs have both energies zero.

[output.txt](output.txt) contains the result data. These checks establish only the
finite range stated above. The universal conjecture requires a general proof;
no Lean proof is supplied here.
