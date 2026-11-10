# Vertex-critical finiteness reversal: exact finite checks

Candidate: Beaton, Cameron and van Omme, *Obstructions to k-colouring H-free
graphs*, arXiv:2610.10737v1, Open Problem 4.3.

Run `python3 check.py`. Python 3.12 or newer and the standard library suffice.
The stored output is `check.out`; successful completion ends in `ALL_CHECKS_OK`.
The program uses integer adjacency bitsets and exhaustive backtracking; it uses
no floating point arithmetic, external graph package, or external program source.

The fidelity check constructs the circulants from Zhou, Jooken, Shan, Goedgebeur
and Huang, *Three-coloring triangle-free graphs without long forbidden paths*,
arXiv:2512.12349v1/v2, Section 4, Theorem 4.1. These are the infinite family cited
in Section 4 of the candidate source. The formula has order 3t+10 and offsets
1, -1, and 5+3j for 0 <= j <= t. For t=3,4,5,6, the program checks triangle
freeness, chromatic number four, every vertex deletion being three-colorable,
and induced-matching number three.

For all 33,868 labeled simple graphs of order at most six, the program checks:

- the Mycielski chromatic increment and induced-matching bound
  nu_ind(M(G)) <= (chi(G)+1) nu_ind(G)+1;
- triangle freeness under lifting, and the triangle-free cutoff
  chi(G) <= 2 nu_ind(G)+1;
- preservation of vertex criticality when chi(G) >= 2;
- the complete multipartite structure of each triangle-containing paw-free
  component, and the clique structure of triangle-containing critical members.

Induced-matching maxima on the base graphs are independently cross-checked by
enumerating vertex subsets whose induced graph has every degree equal to one.
The vertex-recursion solver is used for the lifted graphs. The K1 boundary is
checked explicitly: Mycielski lifting does not preserve its vertex criticality.
The proposed settlement only applies the lifting lemma at chromatic levels
at least four.

For each q=5,...,9, the program also checks the connected complete multipartite
graphs with part sizes (t,1,...,1), for t=1,...,5: chromatic number q,
paw freeness, induced-matching number one, and criticality exactly at t=1.

These are finite checks of the ingredients. They prove neither the infinite
seed theorem nor the quantified open problem, and they do not compute the
existential maximum level m. The proposed positive settlement requires a
general argument and a Lean proof of the published seed family, without a
literature axiom. No Lean proof of the claim or its negation is supplied here.

All implementation code in this directory is original and follows the
repository's license. The mathematical source formulas are attributed above.
