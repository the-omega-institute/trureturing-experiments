# Exact checks for generic HMM identifiability hardness

Source: Markel Zubia and Nils Jansen, *On the Computational Complexity of Hidden Markov Model Identification*, arXiv:2610.09104v1, Problem 36 and Section 8.

Run `python3 check.py`. Python 3.14.8; SymPy 1.14.0. The saved output is `check.out`. All arithmetic is integer, rational, or exact symbolic algebra.

The check reproduces the source's equality between its hidden-path probability formula (Definition 3) and forward procedure (Remark 7) on 7,533 rational cases. It checks equal-emission factorization through word length 5, the symbolic two-marginal inverse for the identifiable family, both selector polynomial maps, and positivity, nonsingularity, distinct-entry rivals and separation on the stated rational grid.

Put p=1/2+u/8 and q=1/2+v/8 for 0<u,v<1. The identifiable family has initial vector (p,1-p), transition rows (q,1-q) and (1/4,3/4), and emission rows (1/4,3/4) and (3/4,1/4). The ambiguous family replaces p by 1/2, q by 1/2+u/8, and both emission rows by (1/4,3/4). Its rivals change u and therefore change a transition entry.

The first two output marginals of the identifiable family are 3/4-p/2 and 5/8+p/8-pq/2. They determine p and q because p>1/2. Its first marginal differs from the ambiguous family's 1/4.

The reductions attach feasibility coordinates to one branch of a tagged union. One uses an identifiable baseline curve and an activated ambiguous component; the other uses an ambiguous baseline curve and an activated identifiable component with a free open square. The dimension comparison is a general mathematical argument about coordinate projections, not a dimension inferred from finite samples.

The output supports a preregistration outline. It does not prove polynomial-time many-one hardness, universal equivalence for all word lengths, semialgebraic dimension theorems, or Lean kernel correctness. The fixed-size finite checks are helpers for the general argument.
