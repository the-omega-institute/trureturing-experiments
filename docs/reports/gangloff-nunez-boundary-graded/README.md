# Residual boundary grading

`check.py` is an original Python 3 program using only integer bit sets and the standard library. Run `python3 check.py`; its deterministic output is stored in `output.txt`. The final line is `ALL_CHECKS_OK` on success.

The source is Silvère Gangloff and Alonso Núñez, *Residual stratification and the Cantor-Bendixson structures of dual algebraic coframes*, [arXiv:2605.04851v3](https://arxiv.org/abs/2605.04851v3), Question 4.20. The source revision is 29 June 2026. Definitions 2.2–2.4, 3.7–3.10, 4.14 and Notation 4.18 determine the tested objects.

The program enumerates all 243 labelled partial orders on zero through four points and constructs their ideal lattices. It calculates co-Heyting subtraction from its defining meet over all admissible ideals. On 1,789 choices of an ideal in these lattices, it checks the source's Lemma 4.2 (the subtraction identity) and Proposition 4.19 (strict monotonicity of the stratum rank). These are fidelity checks of stated source results.

For each integer k from 2 through 10 it checks the poset with b<a and b<c₁<⋯<cₖ. The ideal A={b,a} and the ideal {b} form a cover in the reverse boundary order, with prescribed ranks 0 and k. The k=2 example has seven ideals, derivative sequence P,{b,c₁},{b},∅,∅ and strata {{b,a},{b,c₁,c₂}}, {{b,c₁}}, {{b}}. The output records all reverse boundary covers. The maximal-ideal deletion and principal-ideal subtraction formulas are checked independently against the definition-based calculations.

On the seven-element lattice the program also checks the coframe identity on every subfamily and checks dual compactness against every filtered subfamily, including the empty family as admitted by the source's definition. Precisely the six proper elements are dually compact; top is the meet of the empty filtered family. The output reports 114 filtered subfamilies out of 128 subfamilies.

These are exact finite computations, not Lean proofs. The unbounded family and transfinite stabilization require the mathematical argument in the preregistration; finite checks alone do not prove their universal forms. No third-party implementation or source prose is copied into the program.
