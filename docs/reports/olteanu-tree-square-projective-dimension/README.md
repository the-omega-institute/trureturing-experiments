# Tree-square projective-dimension checks

Run `python3 check.py` with Python 3; only the standard library is required. The exact output is in `output.txt`.

The source is Anda Olteanu, *Edge ideals of squares of trees*, arXiv:2003.10187v2, Question 5.2; Osaka Journal of Mathematics 59 (2022), 369–386, Question 6.2, DOI 10.18910/87483. The question compares the projective dimensions of the edge ideals I(T) and I(T²), where T² is a graph distance square, rather than an ordinary ideal square.

The program computes projective dimensions of the quotient S/I(G) by Hochster's formula. It enumerates every induced independence complex and computes its reduced simplicial boundary ranks by exact Gaussian elimination over the rationals, F₂ and F₃. For nonzero ideals it subtracts one to obtain pd I(G). For zero ideals it uses a separate bottom value, without subtracting one. Distances come from breadth-first search; unreachable vertices remain unreachable.

Fidelity: the quotient projective dimension of the square of the n-vertex path is n − ceil(n/5), as given by Proposition 4.3 and Corollary 4.4 of the preprint. This is reproduced for 2 ≤ n ≤ 8. The Morey formula recalled in its final section, pd S/I(P_n) = n − ceil(n/3), is checked on the same sizes.

Settlement screen: the comparison holds for all 1,100 labelled simple graphs with 0 ≤ n ≤ 5 and all 1,442 labelled trees with 1 ≤ n ≤ 6, for distance powers 2 ≤ r ≤ max(2,n−1), over each of the three fields. For graphs without isolated vertices, the program also checks the proposed domination and codimension inequality chain. The Prüfer-word enumeration includes every labelled tree exactly once at each size. Empty graphs, isolated vertices and disconnected graphs are included in the simple-graph screen.

All arithmetic used as evidence is integer or rational, with modular inverses in finite fields. These are finite tests only: they prove neither the universal claim nor field independence. The planned positive settlement uses a general argument, rather than finite enumeration.

The script is original research-check code. Its mathematical inputs are the source's formulas, Hochster's formula, and Dao–Schweig, *Projective dimension, graph domination parameters, and independence complex homology*, arXiv:1110.2841, Theorem 4.4 and the packing observation in the proof of Corollary 6.12.
