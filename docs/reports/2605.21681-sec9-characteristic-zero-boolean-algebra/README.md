# Boolean incidence kernels for the characteristic-zero finite-length question

Source: Jingjie Yang, Mikołaj Bojańczyk and Bartek Klin, *The Finite Length Property of the Rado Graph and Friends*, arXiv:2605.21681v2, Section 9; LICS 2026, article 82, DOI [10.4230/LIPIcs.LICS.2026.82](https://doi.org/10.4230/LIPIcs.LICS.2026.82). The source asks whether every oligomorphic structure has the finite length property over characteristic-zero fields and explicitly names the countable atomless Boolean algebra as an unresolved example.

Run `python3 check.py` with Python 3.10 or later; only the standard library is required. The recorded output is in [output.txt](output.txt). Arithmetic is integral, including fraction-free row elimination computing ranks over the rationals. No floating point is used.

The fidelity check reproduces the five ordered-pair orbit types of vector atoms over the two-element field (Section 2, Example 2.5), using all 168 elements of GL(3,2). The orbit sizes are 1, 7, 7, 7 and 42. Every infinite-dimensional pair type is already represented in dimension three; every finite linear automorphism extends to the infinite-dimensional space.

For the Boolean construction, an algebra with n atoms is represented by bit masks, excluding the empty and full masks from the basis. The row indexed by a tuple tests whether all tuple atoms lie in a basis subset. A row depends only on the nonempty set of distinct tuple entries; repetition realizes every such set of size at most d. The program computes the ranks of these integer incidence matrices over the rationals. Restriction to the first k binary coordinates gives n = 2^k atoms and an algebra of size 2^(2^k).

| k | n | Dimension of the free space | Dimensions of kernels K_d for d = 1, ..., n-1 |
|---|---|---|---|
| 1 | 2 | 2 | 0 |
| 2 | 4 | 14 | 10, 4, 0 |
| 3 | 8 | 254 | 246, 218, 162, 92, 36, 8, 0 |

The program checks closure of the row family under adjacent atom transpositions, all refinement identities for surjective maps between atom sets with at most four fine atoms at d = 1, 2, 3 (50,648 basis/tuple equalities), and surjectivity on tuple sets, which makes pullback injective.

For d = 1, ..., 5 it exhausts all ordered d-tuples in a partition with d+3 atoms. The alternating vector on the 2^(d+1) proper nonempty subsets containing a fixed atom c has every d-tuple test zero and has test one at the (d+1)-tuple of distinct variable atoms. Integer identities also verify the corresponding readings modulo 2, 3, 5 and 7.

These are finite checks supporting a proposed argument, not a proof of the infinite statement. Countability, atomlessness, extension of finite Boolean isomorphisms, invariance under the infinite automorphism group, compatible directed unions and arbitrary chain length remain mathematical obligations. The intended mechanism applies to every field and refutes the characteristic-zero assertion already over the rationals; it makes no assertion about ascending chains or homogeneous structures in a finite relational vocabulary.
