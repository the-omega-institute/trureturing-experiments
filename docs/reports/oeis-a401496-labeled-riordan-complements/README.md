# Labeled Riordan graphs with labeled Riordan complements

Run `python3 check.py` with Python 3. Only the standard library is required.
`output.txt` is the captured output. All arithmetic is exact.

The definition is the first comment of [OEIS A401493](https://oeis.org/A401493):
for source labels I > J, adjacency is the coefficient of t^(I-2) in g f^(J-1)
over GF(2), with f(0)=0. For each n from 1 through 10 the program enumerates
all truncated pairs, constructs and deduplicates adjacency matrices, and tests
complements with the labels fixed. As fidelity checks it reproduces the total
labeled Riordan count (4^(n-1)+2)/3, all available A401496 terms in this range,
and Pascal adjacency by independently computing binomial parity.

A separate construction enumerates the proposed classification using
r delta^2=0 modulo t^(n-1), and compares the resulting graph sets exactly to
the exhaustive ones. It checks the injective quotient-class counts at every r;
it does not count duplicate presentations as distinct graphs. The finite sum,
recurrence, differences, A113836 identity, parity formulas and coefficients
of the rational generating function agree through index 200. These checks are
finite diagnostics and supply no all-order proof or Lean proof probe.

Sources: Arnauld Mesinga Mwafise, [OEIS A401496](https://oeis.org/A401496),
revision 9; A401493, revision 10; Artur Jasinski, [A113836](https://oeis.org/A113836).
Cheon–Jung–Kitaev–Mojallal, *Riordan graphs I: Structural properties*,
[arXiv:1710.04604v2](https://arxiv.org/abs/1710.04604v2), is context.
Mwafise's [Zenodo 23676629](https://doi.org/10.5281/zenodo.23676629) supplies
finite counts through n=14 and explicitly does not supply mathematical proofs.
The program is independently written from the coefficient definition; no code
from that software is incorporated.

The parity FORMULA line in OEIS HTML, JSON and revision history literally prints
`108^(k-1)` and `38^k`. The intended products `10*8^(k-1)` and `3*8^k` agree
with the other four conjectures and all 14 data terms. Reading the concatenations
as integer bases instead gives a(3)=6 and a(4)=16, whereas direct counts are 4
and 12. The program reports this source discrepancy explicitly. Its successful
identity checks use the disclosed multiplication reading, and do not propose
settling the typographical error by refutation.
