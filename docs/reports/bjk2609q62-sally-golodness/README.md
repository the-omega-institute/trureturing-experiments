# Sally two-omission family: exact finite checks

Source: Manav Batavia, Omkar Javadekar and Manohar Kumar, *The Betti numbers and
Golodness of numerical semigroup rings of Sally type*,
[arXiv:2609.26752v1](https://arxiv.org/abs/2609.26752v1), Theorem A(b,c) and
Question 6.2(a)–(c).

Run `python3 check.py`. The program uses only the Python standard library,
integer arithmetic and exact elimination over the prime fields F₂ and F₃.
Its output is in `output.txt`; successful completion ends in `ALL_CHECKS_OK`.

The semigroup generators determine the Apéry representatives by integer
enumeration. Their sums determine the multiplication of A = k[S]/(tᵐ),
independently of the proposed multiplication table. The program checks that table
for all 66 pairs 5 ≤ m ≤ 15, 4 ≤ n < m.

For (m,n) = (5,4), (6,5), (7,6) over both fields:

- Koszul homology on the degree-one generators of A reproduces the source's
  defining-ideal Betti polynomials, with coefficient lists (1,3,2), (1,6,8,3)
  and (1,9,17,12,3), respectively.
- The reduced normalized bar complex computes β₀,…,β₄ of k over A by independent
  exact matrix ranks. The program also checks d² = 0 over the integers.
- The standard regular-linear-element identity P_R = (1+z)P_A transfers those
  numbers to k over R. Thus the reported R lists use this mathematical identity;
  they are not ranks of a directly truncated bar complex over the infinite ring R.
  Here tᵐ is regular because R is a polynomial subring and is outside the square
  of the augmentation ideal because every positive semigroup element is ≥ m.

The series check tests the proposed first Serre defect and the full difference
identity through degree 30 on the same 66 parameter pairs. A separate quotient
rank calculation uses just the two defining quadratic-dual relations, and
compares its dimensions in degrees 2–4 over both fields with the proposed
normal-word count. The latter is also checked through degree 6 on five pairs,
including n=6 and n=m−1.

These finite checks support the proposed formulas. They do not prove them for
all parameters, all homological degrees or all characteristics, and do not
establish a Lean settlement. The characteristic-free minimal resolution, the
regular-element transfer and the identifications with Tor remain proof
obligations of that settlement.
