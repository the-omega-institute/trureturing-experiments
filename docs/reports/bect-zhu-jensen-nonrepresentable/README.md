# Bect–Zhu Jensen functionals without a loss representation

Source: Julien Bect and Xujia Zhu, *Uncertainty functionals revisited:
Concavity and Jensen’s inequality*, arXiv:2605.18100v1, Open Problem 2.
The fidelity target is Section 2.2.1, Proposition 2, equations (5)–(6),
specialized to the Gini–Simpson entropy h(t)=t−t².

Run `python3 check.py` with Python 3 and SymPy 1.14.0. All evidence uses
integers, fractions.Fraction, and exact symbolic algebra. The captured standard
output is check-output.txt; success requires exit 0 and final line ALL_CHECKS_OK.

The candidate is H(p)=min(1,liminf 2ⁿpₙ) on all probabilities on the discrete
natural numbers. For 0<c<1, p⁽ᶜ⁾ has p₀=1−c and pₙ=c2⁻ⁿ for n≥1. The perturbation
q⁽ᶜ,ᴺ⁾ has q₀=(1−c)(1−2⁻ᴺ), qₙ=c2⁻ⁿ for 1≤n≤N, and qₙ=2⁻ⁿ for n>N.
It has H(q)=1, while H(p)=c. Profiles in the program store an exact geometric
infinite tail, including its entire mass; no tail is silently discarded.

The checks reproduce the source’s finite-mixture Jensen-gap identity, test
normalization, full support, domination, clipping and Dirac values, and check
an exact risk-difference formula for several infinite integrable loss families.
Symbolic identities verify normalization and geometric-tail moment formulas.
The mathematical argument for arbitrary nonnegative losses uses the vanishing
tail of a convergent nonnegative series, rather than enumeration of losses.

These results support a preregistration. They do not prove measurability, Jensen
for arbitrary random probability measures, or nonrepresentation for every loss
family. Those assertions require the general argument and subsequent Lean proof.

Program author: codex-cli research flight d3a4595723089036ff47ce60, using the
search-seat candidate supplied in the research brief. Program license: MIT.
