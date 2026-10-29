# Balanced complete multipartite general-position polynomials

Run `python3 check.py` with Python 3.10 or later. The program uses only standard-library exact integer and rational arithmetic and prints `ALL_CHECKS_OK` on success. The captured output is in `output.txt`.

Source: Bilal Ahmad Rather, *Explicit Formulas and Unimodality Phenomena for General Position Polynomials*, arXiv:2603.06930v3, Sections 3–4 and Section 8 item 2. Here t is the number of parts and n is their common size, called r in the source.

The program reproduces the source's coefficient lists for K_(8,8) and K_(5,5), and checks its small-part log-concavity theorem for 2 ≤ t ≤ 100 and 1 ≤ n ≤ 4. An independent enumeration uses BFS distances and forbidden geodesic triples on every vertex subset of 14 graphs, totaling 108716 subsets, to check the counting formula and polynomial identity.

For 2 ≤ t ≤ 15 and 1 ≤ n ≤ 400, all 5600 parameter pairs agree with the proposed classification:

- Log-concavity holds exactly when n ≤ t+1 or (t,n)=(2,4).
- Unimodality holds exactly when n ≤ 2t+3 or t C(n,t)(n−2t−1) ≥ (t+1)n^t.

The same rectangle has no log-concavity failure outside index t+1. Exact junction comparisons check the proposed recovery-threshold boundaries and the strict-valley family n=2t+4 for 2 ≤ t ≤ 15. Threshold searches use the proposed monotonicity of the junction expression; they check the two neighboring integers directly and do not prove that monotonicity.

These finite checks corroborate a proposed uniform mathematical argument. They do not prove the classification for unbounded t or n and contain no Lean proof.
