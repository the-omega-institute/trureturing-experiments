# Threshold-three torus winding

The reference rule in Muñoz–Nuño, *Symmetry-resolved periods in a cyclic-dominance automaton*, [arXiv:2610.37990v1](https://arxiv.org/abs/2610.37990v1), advances a color when at least three of its eight Moore neighbors have the successor color. Section V.4 leaves conservation of torus winding for future work.

Run `python3 check.py`. The program uses integers and exact rational comparisons. Updates enumerate physical neighboring sites independently of the residue argument. Its stored standard output is `check-output.txt`; the final line is `ALL_CHECKS_OK`.

Fidelity: all 220 positive three-band width triples for square sizes 3 through 12 translate by one site, preserve populations, and change exactly three columns, as described in Section V.4. The equal-width size-six example has exact state period six. Exact comparisons also reproduce equivalence of the source parameters α = 1/3 and α = 3/8.

Candidate refutation: for L = 4m and k = (j − i) modulo 4, use the profile (2, 1, 0, 2). For m = 1 through 16 and thresholds three and four, one update yields the fixed profile (2, 2, 0, 2). All rows change winding from −m to zero, all columns from +m to zero, and all elementary-square winding sums vanish before and after. Every initial field contains all three colors. The output includes the complete size-eight state, predator counts, updated state, loop windings, and square sums. All 81 four-residue color profiles are also checked for zero square sums.

These are finite exact computations. The infinite-family argument is part of the preregistration draft; no Lean proof or kernel verification of the settlement is claimed by these files.
