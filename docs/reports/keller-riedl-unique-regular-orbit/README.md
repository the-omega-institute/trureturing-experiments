# Keller–Riedl Question 9.3: exact checks

Run `python3 check.py` with Python 3. The program uses only standard-library integer arithmetic; its deterministic output is in [output.txt](output.txt). A successful run exits zero and ends with `ALL_CHECKS_OK`.

The reference is Keller–Riedl, *Linear 2-Group Actions With Exactly One Regular Orbit*, arXiv:2610.27536v1, Sections 2, 4 and 5. Fidelity checks reproduce the Singer-type element, the two local regular-orbit stabilizers (Lemma 5.3), all binary-tree subset orbit shapes through level five (Theorem B/4.5), and the actual permutation groups R₃ and R₄ (Corollary 5.11).

The proposed negative answer replaces each eight-coordinate local group Σ by an index-two parity kernel in its extension by a simultaneous Frobenius/permutation involution. The program checks its normalizer generators, exact composition and transporter formulas, local stabilizers, transitivity, and all 80,640 permutation/global-Frobenius candidates in the local diagonal centralizer. The resulting local center consists of the two prime-field scalars. The order identity for the lifted family is checked at j = 3,…,64.

These are finite checks of ingredients, including some exhaustive local checks and some sampled checks. They do not prove closure, regularity, or uniqueness for every group element, do not enumerate the first dimension-128 example, and do not prove the infinite family or any Lean theorem. The symbolic argument must establish those conclusions. The proposed family has q = 3, m ≥ 8, order 2^(99·2^(m−7)−1), and center of order two; every source T(R) has a central element of order eight. The q ≥ 5 classification remains outside this construction.
