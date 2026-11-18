# Exact checks for Nies Problem 2.1

The target is André Nies, Problem 2.1 in *Open Problems in Mathematical Logic*, arXiv:2608.26628v2: finiteness of the continuous outer automorphism group of the automorphism group of every countably infinite omega-categorical structure in a finite first-order language.

Run `python3 docs/reports/nies2608p21/check.py` from the repository root. The program uses only Python's standard library and exact integers and rational fractions. `output.txt` is its complete successful output.

The fidelity check reproduces the source's displayed order-reversal representative on rational increasing affine maps: conjugation by negation sends `(a,b)` to `(a,-b)`, is an involutive homomorphism, and changes translation by +1 to translation by -1. It checks 253 maps and all 64,009 composition pairs. It also enumerates the two possible actions on the pair-order types equality, less-than and greater-than. This verifies the displayed representative calculation on an exact rational subfamily; it does not independently classify all continuous automorphisms of Aut(Q,<).

The settlement checks exhaust every permutation of B×C2 for base sizes 1 through 4 and verify the synchronized-cover automorphism description. They check every shear, every input pair and every composition of shears for C2×C2^d with 0≤d≤6. They also check C2×(C2^d×S3), 0≤d≤3, enumerating all inner automorphisms and testing equality of outer classes. Nonempty characters have explicit coordinate witnesses. No finite test establishes the unbounded family or the existence, omega-categoricity or topology of the infinite finite-language coding. Those are obligations of the proposed general argument.

The mathematical ingredients are Bodirsky–Evans–Kompatscher–Pinsker, arXiv:1510.00356v2, Propositions 3.4 and 3.18 and the paragraph immediately following Proposition 3.18. The proposed combination transfers independent continuous C2 characters through dense restriction and realizes them as central shears of the synchronized cover. The program is independently written and contains no third-party program or Lean source.
