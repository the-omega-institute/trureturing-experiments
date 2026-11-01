# Directed boundary cycle-cover checks

This self-contained Python 3 program uses integers and `fractions.Fraction` only.
Run from this directory:

```sh
python3 check.py
```

`output.txt` contains the complete successful output and ends in `ALL_CHECKS_OK`.

The source is Mohamed Ali Belabbas, *Structural stability of systems and cycle
covers in random graphs*, [arXiv:2610.01607v1](https://arxiv.org/abs/2610.01607v1),
Definitions 5 and 7–9, Lemma 7, Proposition 1, and Section 9's first open question.
The graph model allows loops and independently samples the two orientations.
Property N requires a cycle-covered induced subgraph of every cardinality;
those subgraphs need not form a nested sequence.

The fidelity check independently compares the type-level Hall inequalities with
bitmask perfect-matching search in the complete directed blow-up, for every
strongly connected skeleton with two or three types and every type-count vector
of total size 1 through 6: 148 skeletons and 12,060 vectors.

The proposed boundary family has arcs 1→1, i→i+1 for i<q, and j→1 for j≥2.
The program enumerates its simple cycles, checks the exact inverse barycentric
formula for its ordered-coordinate polytope, its affine dimension, its active
Hall inequalities at the uniform concentration vector, and the categorical
covariance. It verifies one accepted ordering chamber out of q! for q=2 through
6. These rational chamber fractions are 1/2, 1/6, 1/24, 1/120 and 1/720; they
are not Gaussian probabilities estimated by simulation.

An independent cycle-packing recursion checks the whole-graph and all-sizes
properties for every count vector with q=2 through 5 and total size 1 through 10
(4,352 vectors), when all supported edge probabilities equal one. Exact
multinomial weights verify normalization and the tie-free chamber identity.
For q=2 they also reproduce the symmetry/atom formula used in the source's
Proposition 3: Pr(N1≥N2)=1/2+Pr(N1=N2)/2. The no-loop singleton obstruction
and a looped dimension-deficient triangle check the two proposed zero-limit
mechanisms.

These are finite consistency checks. They do not establish the unbounded
asymptotic limit, a uniform probabilistic matching estimate, a categorical CLT,
or the claim for arbitrary supported edge probabilities. The planned positive
settlement uses a general proof, not finite enumeration.

The check program is original code and contains no copied third-party code.
