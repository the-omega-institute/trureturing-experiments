# Hindi–Schutzer endpoint nonlocality: exact checks

Source: Felipe de Mattos Chafik Hindi and Waldeck Schutzer, *Atomic
Decompositions of Lie Characters and the Dominant Weight Poset*,
[arXiv:2609.08024v1](https://arxiv.org/abs/2609.08024v1), Question 6.12.
Bourbaki numbering is used throughout; endpoint restriction preserves Dynkin
labels as in source equation (22).

Run `python3 check.py`. The program uses only the Python standard library,
with integer and `fractions.Fraction` arithmetic. Its saved output is
`output.txt`; success ends with `ALL_CHECKS_OK`.

The program computes weight multiplicities by Freudenthal's recurrence with
Weyl normalization, enumerates the dominant interval, and evaluates its
Möbius function directly. It independently inverts the triangular character
expansion, enumerates all atoms and all their joins, and checks the crosscut
sum. For ranks 4–6 and s = 0, 1, it also evaluates the full signed Kostant
Weyl sum at every weight of the interval and checks agreement with
Freudenthal, including vanishing outside the A-type parabolic fixing the
first orthogonal coordinate. A third partition calculation removes the
parabolic roots, and a multigraph enumeration verifies the resulting
partitions with at most two edges.

Fidelity checks reproduce the source's D4 coefficients
`a(0100,0022) = -2` and `a(1100,1022) = -1`.
For every 4 ≤ n ≤ 8 and 0 ≤ s ≤ 3, the program checks
μ = ω_(n−2) + sω_n, λ = 2ω_1 + (s+2)ω_n,
and the proposed formula `2a(μ,λ) = n²−9n+16`.
It checks the three prescribed atoms, all join multiplicities, and the
sharper comparison λ′ = (s+2)ω_n, whose ambient coefficient is zero.
For n ≥ 5, all six identifications of the unique radius-one D4 diagram give
zero for the projected pair. Thus D5 gives −2 and D6 gives −1 despite local
coefficient zero; λ and λ′ have identical local endpoint labels.

These finite computations support a proposed refutation and the family
formula. They are not a Lean proof or a proof for unbounded n and s.
The mathematical definitions and algorithms were implemented from the
paper and the stated recurrence formulas; no external program source was
copied. This report's code is offered under the repository's license.
