# Orbital entropy partition: exact preregistration checks

The target is Conjecture 2 of Damiano Aliverti-Piuri, Julia Liebert and Christian
Schilling, *Maximal Orbital Entanglement from Non-Interacting Fermions*,
arXiv:2610.32974v1. The claim concerns subset-average patterns of the binary
spectral-entropy maximizer for diagonal block spectra. It imposes no trace-one
normalization.

Run `python3 check.py`; the expected output is stored in `output.txt`. Python 3
and SymPy are required. All evidence uses integers, `fractions.Fraction`, or exact
SymPy algebra. No floating-point tolerance, sampling, or numerical optimizer is
used.

- Source fidelity: reproduce the characteristic polynomial and block spectrum
  of the source's displayed 3-by-3 optimizing matrix for the partition (2,1),
  its Definition 4 averaging example, and the binary-entropy endpoint values.
- Exhaustive small input cases: occupation entries are multiples of 1/8 in
  [0,1], in non-increasing order. All admissible partitions for d=2,3 give 375
  input/partition cases. For the finest partitions, exact Fourier unitaries
  attain the Jensen upper bound. For (2,1), the feasible region is the source's
  full Cauchy-interlacing rectangle. Intersect it with all six coordinate-sorting
  chambers and enumerate every vertex with rational linear algebra. The top-k
  sums are affine on each chamber, so their inequalities at its vertices
  certify the whole continuous region, including degenerate rectangles.
  The 865 vertices satisfy the majorization inequalities for the proposed
  maximizer. Strict concavity of binary entropy gives uniqueness. Enumerating
  set partitions independently finds the required averaging witness.
- Endpoint family: all 50 triples with 2 <= d <= 8,
  ceil(d/2) <= m < d, r=d-m, and r <= N <= m are checked for
  lambda=(m,1,...,1), n=(1^N,0^(d-N)). Exact rank-N projections attain the
  bound 2r log(2). The large block must have at least N-r eigenvalues 1 and m-N
  eigenvalues 0, by intersection/rank bounds, leaving at most 2r nontrivial
  eigenvalues in the full block spectrum. Each contributes at most log(2).
  The program verifies characteristic polynomials, idempotence, traces,
  commutation, entropy equality, and the disjoint averaging partition.
- Five exact matrix examples verify the sign in the proposed unitary
  first-variation commutator calculation.

These checks certify the stated finite cases using standard interlacing,
majorization, concavity and rank facts. They do not prove the all-dimensional
conjecture, the analytic spectral derivative lemma, uniqueness in general, or
Conjecture 1. No Lean proof of the claim or its negation is part of this report.
