# KMV Question 12: direct-sum heredity

Exact evidence for the proposed negative answer to the direct-sum clause of
Kaygorodov–Martín González–Viruel, *Evolution algebras: a survey of results and
a list of open problems*, arXiv:2610.22248v1, Question 12. The characterization,
graph, tensor-product and ideal clauses are outside this report.

Run `python3 check.py` with Python 3 and SymPy (checked with SymPy 1.14.0).
`output.txt` is its complete output. Successful execution exits 0 and ends with
`ALL_CHECKS_OK`. Arithmetic uses integers, `fractions.Fraction` and exact SymPy
expressions; floating point is not used.

- Source fidelity: Section 1.1's first example, the function algebra K^n has
  2^n idempotents, is reproduced for n=1,...,5 over F_3 and F_5.
- Family: m=4,...,9 and n=3,...,8 over Q, F_3, F_5 and F_7; the factors have
  one-dimensional product spans and zero annihilators. The four proposed
  generators are independent and all 16 ordered products match the table.
- All subspaces of A_4 and B_3 over F_3 and F_5 are enumerated once using RREF.
  For every closed subspace, the program constructs and verifies an intrinsic
  natural basis. There are respectively 50, 10, 110 and 14 closed subspaces.
- The quotient multiplication formula is checked symbolically. All ordered
  pairs in the two-dimensional quotient over F_3, F_5, F_7, F_11 and F_13 are
  checked: no independent pair has zero product.

These are exact finite and symbolic checks. The all-field, unbounded-dimensional
family and the absence of an arbitrary natural basis require the mathematical
argument in the preregistration and a later Lean proof. This report contains no
Lean proof or claim of kernel verification. Literature priority is not tested by
the program.
