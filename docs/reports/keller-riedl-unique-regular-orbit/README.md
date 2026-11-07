# Keller–Riedl Question 9.3: exact checks

Run `python3 check.py` with Python 3. The program uses only standard-library integer arithmetic; its deterministic output is in [output.txt](output.txt). It exits zero and ends with `ALL_CHECKS_OK` when the fidelity checks and the failure diagnosis below are reproduced. This does **not** mean that the expected refutation passed.

The reference is Keller–Riedl, *Linear 2-Group Actions With Exactly One Regular Orbit*, arXiv:2610.27536v1, Sections 2, 4 and 5. Fidelity checks reproduce the Singer-type element and the two local regular-orbit stabilizers (Lemma 5.3), every binary-tree subset orbit shape through level five (Theorem B/4.5), and the actual permutation groups R₃ and R₄ (Corollary 5.11).

The candidate construction replaces the eight-coordinate group Σ by the parity kernel S in its extension by θ, where θ applies Frobenius in all coordinates and permutes labels by multiplication by three. Exact normalization, composition, transporter and stabilizer checks support those local ingredients.

The proposed center obstruction fails. The exact center calculation screens all 80,640 local permutation/global-Frobenius pairs and all 64 remaining even/odd diagonal candidates. It finds eight central elements, with diagonal exponent pairs (a,3a mod 8), rather than two. In particular the diagonal transformation with ζ on even coordinates and ζ³ on odd coordinates belongs to S and is central of order eight. Commuting with a semilinear permutation relates its coefficients by Frobenius; transitivity of the permutation projection does not force them to be constant. Repeating this diagonal transformation in every block also gives a central element of order eight in every proposed U_m.

Thus the asserted centers Z(S)=C₂ and Z(U_m)=C₂, and the proposed nonconjugacy argument based on them, are false. The calculation does not establish whether these groups are conjugate to any T(R), and does not settle Question 9.3. A proposed preregistration relying on this argument must abstain. No Lean proof is provided.
