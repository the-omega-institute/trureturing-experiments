# Star cluster component counts by codimension

Exact checks for Eleonore Faber and Bernd Schober, *Singularities of star cluster algebras*, arXiv:2609.01501v1, Proposition 4.3, Corollary 4.4 and Question 4.5.

Run `python3 check.py` with Python 3.12 or newer and Singular 4.4.1 available on PATH. The stored result is `output.txt`. All evidence uses integer or rational arithmetic and polynomial ideals over the rationals and the field with two elements.

The check covers all 117 ordered tuples with 2 ≤ n ≤ 4 and 1 ≤ q_a ≤ 3, comprising 306 stage distributions. Each tuple is explicitly identified with its sorted representative by a block-coordinate permutation. Singular computes the minimal prime ideals of the source's defining zero-block/difference ideals for all 31 sorted representatives, all permitted stages and both coefficient fields: 142 stage decompositions. Every returned prime is matched by equality of Gröbner ideals to exactly one predicted component, and its ambient codimension is computed independently by Singular. Intersecting the source stage ideals also verifies that the literal union has the last stage's ideal.

Separate finite enumeration checks the generating polynomial, multiplicity-based counts, exact support, noncontainment within each stage, containment in the last stage, the pair formula for the literal union, and the source's unrefined counts. The tuple (1,1,4,4) illustrates support gaps; its literal union has counts 1, 8 and 4 at codimensions 6, 12 and 18.

These are finite consistency checks. They do not prove the universal assertion over every field or constitute Lean verification. No implementation of a third-party proof is copied.
