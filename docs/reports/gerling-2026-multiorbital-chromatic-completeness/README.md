# Gerling multiorbital chromatic completeness

`check.py` uses Python 3 and SymPy for exact integer enumeration and rational polynomial algebra. Run:

```sh
python3 check.py
```

`output.txt` contains its output. Success requires exit 0 and the final line `ALL_CHECKS_OK`.

The fidelity check reproduces the reflected three-vertex path and the edge-plus-isolate example in Melanie Gerling, *The Multiorbital Bivariate Chromatic Polynomial: A Subgroup-Lattice Refinement*, arXiv:2609.01866v4, subsection 3.4: the orbital polynomials agree and the multiorbital difference is `y(1-x)`.

The proposed completeness counterexample starts with the two five-vertex edge sets in `check.py`. The program checks their common polynomial and both palette-conditioning identities. Vertex-partition enumeration independently computes the polynomials of their joins with complete graphs for `0 <= r <= 5`; for `2 <= r <= 5` it checks the quotient under the clique transposition, equality of the multiorbital polynomials, connectedness, faithful C2 actions, and different universal-vertex counts. Full inclusion-exclusion checks all 131072 edge subsets of each seven-vertex graph. Direct enumeration checks every coloring and every transposition-fixed coloring for all 15 palettes `0 <= y <= x <= 4`, on both base graphs and both seven-vertex graphs.

The finite calculations corroborate the proposed uniform family. They do not prove the universal family statement or replace its counting-polynomial bridge in Lean. The complete mathematical outline belongs to the preregistration draft. No Lean proof is included here.
