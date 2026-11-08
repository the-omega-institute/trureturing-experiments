# Sparse output controllability: exact finite checks

The proposed characterization uses the early delays 0 through d−1 modulo the output
span of the stable reachable tail. It predicts sparse output controllability exactly
when every cut J satisfies rank([C A^d W | C A^j B, j∈J]) + s(d−|J|) ≥ p.
The check is evidence on finite rational systems; the universal statement requires
a separate mathematical proof and Lean verification.

Run `python3 check.py`. Python 3.12 or newer suffices; no external package is used.
All arithmetic uses integers and `fractions.Fraction`.

- Source fidelity: Joseph, arXiv:2005.12746v2, Examples 1 and 2. Both ordinary
  output ranks, maximal rank drops (2), and maximal prefix averages (1) agree.
  Exhaustive full-horizon schedules give ranks 2<3 and 2=2, respectively.
  Example 2 also has rank 2 using only the second actuator, as stated in the source.
- Seed 377: 120 random integer systems, d=1,…,4, m=1,…,3, with dense, singular,
  and strictly upper triangular A. The cut minimum equals the maximum early
  selected-column rank augmented by the tail in every system: 45 feasible and
  75 infeasible. Each feasible case receives a constructed schedule at d(r+1),
  solves an exact rational target from a nonzero initial state, and verifies it
  by the physical recurrence.
- Every schedule at the asserted horizon is enumerated for the 112 random systems
  with at most 4096 maximal support schedules: 8380 schedules in total. Supports
  of size less than s need not be separately enumerated because they extend to
  size s and extension cannot decrease column rank.
- For the 98 random systems with dm≤8, all 4364 ground-set subsets are enumerated.
  The unrestricted intersection cut, compressed delay cut, and maximum common
  independent-set size agree. This checks the cut-compression step independently
  of whole-horizon enumeration.
- Nilpotent-chain timing obstruction: k=2,3,4. Only delay 1 contributes; direct
  two-delay enumeration gives minimum sparsity k. Since A²=0, later delays add
  zero columns. The original necessary integer bound gives only ceil(k/2).

Small source ranks are independently checked by determinant minors. Output is in
`check-output.txt`, ending in `ALL_CHECKS_OK`.

Sources: Geethu Joseph, *Output Controllability of a Linear Dynamical System With
Sparse Controls*, arXiv:2005.12746v2, DOI 10.1109/TCNS.2022.3188484; Alexander
Schrijver, *A Course in Combinatorial Optimization*, §10.5, Theorems 10.8–10.9,
https://homepages.cwi.nl/~lex/files/dict.pdf. Source texts are cited, not distributed.
