# Evolution-algebra automorphism realization checks

Source: I. Kaygorodov, C. Martín González and A. Viruel, *Evolution algebras: a survey of results and a list of open problems*, arXiv:2610.22248v1, Proposition 14 and Question 31. The proposed negative answer concerns full algebraic groups and group schemes over an algebraically closed field of characteristic zero. Abstract-group realization is outside this check.

Run `python3 check.py` with Python 3 and SymPy. `output.txt` is the resulting output. Arithmetic is exact throughout; successful completion ends in `ALL_CHECKS_OK`.

The program reconstructs all Leibniz equations of Proposition 14's algebras for 1 ≤ m ≤ 7 and verifies that their kernels equal the skew-symmetric m-by-m matrices extended by zero. It checks nontrivial rational orthogonal automorphisms, generic dual-number multiplication and determinant identities in dimensions 1–3, and the nondegenerate adjoint Killing form of sl₃.

For the proposed SLₙ obstruction it verifies the two polynomial identities obtained from equality of ranks and dimensions. It also checks 2 ≤ n ≤ 256 and records the exceptional invariant matches n=2 and n=4. Equal dimensions alone have other matches, so the full-family argument requires the rank invariant.

These checks support a preregistration argument. They do not prove the classification of derivation Lie algebras, differentiation of natural group-functor isomorphisms, the general Cartan-rank formulas, or non-realizability for arbitrary evolution algebras. The planned proof reuses the repository's frozen semisimple-derivation classification. No Lean proof of Question 31 or its negation is included.
