# Exact subgroup lower-bound checks

Source: Stefanos Aivazidis, *Subgroup bounds for abelian p-groups with
applications*, [arXiv:2606.03320v1](https://arxiv.org/abs/2606.03320v1),
Section 2 and Question 10.1. The script is independent code implementing the
published enumeration formula; no source text or third-party code is included.

Run `python3 check.py` with Python 3. Only the standard library is required.
The stored output is `output.txt`. All arithmetic uses integers, integer
coefficient polynomials, and `fractions.Fraction`.

For a partition λ, write d_h = λ_h − λ_(h+1), q = 1/p,
η = Σ_i floor((λ'_i)^2/4), and θ = ∏_(h odd)(d_h+1).
The proposed normalized polynomial lower bound is

\[
F_\lambda(q)=(d_1+1)\prod_{a\ge1}
\left(1+d_{2a+1}[a+1]_q+[a]_q\sum_{j=1}^{d_{2a}}q^j\right),
\qquad [m]_q=\sum_{j=0}^{m-1}q^j.
\]

The proposed height factors are β_(2a) = 2c_a/(1+c_a) and
β_(2a+1) = (1+c_a)/2, where c_a = [a+1]_q and a ≥ 1.
They multiply θp^η once per occupied height h ≥ 2.

The program checks:

- Fidelity: the complete polynomial for λ = (5,5,3,1) printed in the source.
- The leading coefficient θ and degree η, and the coefficientwise inequality
  p^(−η)s(A_λ) ≥ F_λ(q), for every partition of size at most 12 and the
  source example. It also checks that F has a positive next coefficient
  for every noncyclic partition in that range.
- Validity and distinctness of all selected paired-height subgroup types
  in those cases, including the empty partition and empty height blocks.
- Both proposed inequalities at p = 2,3,5,7 using exact rational arithmetic.
- The Delsarte counts against a separate enumeration of actual subgroups,
  using only addition in products of residue groups, for all partitions of
  size at most 5 at p=2, at most 3 at p=3, and at most 2 at p=5.

The checks cover 273 partitions, 1,092 rational cases, 2,875 selected types,
and 30 independent subgroup enumerations. They pass with final line
`ALL_CHECKS_OK`. These are finite checks of a proposed general argument;
they do not establish the assertion over every prime and partition.
