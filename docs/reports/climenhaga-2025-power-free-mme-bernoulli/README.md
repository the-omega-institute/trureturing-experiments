# Power-free MME Bernoulli proposal: exact checks

This report checks finite helper claims for a proposed affirmative answer to
Vaughn Climenhaga, *Maximizing entropy for power-free languages*,
arXiv:2507.18779v2, Question 1.4 in the published version
(DOI 10.1017/etds.2026.10359).

Run `python3 check.py` with Python 3 and SymPy. The saved output is `output.txt`.
All evidence uses integers, rational numbers, or symbolic polynomials.

The source-fidelity checks reproduce the source's fourth-power-end count
`|G_4| = d^4 - d` and its binary length-eight lower bound of 188. The exact
length-eight count is 194. For the proposed product proof, the program checks
the binary length-sixteen count 49100 and union lower bound 46752, the strict
entropy-gap inequalities after exponentiation, the composition bound, local
gluing on the reported finite ranges, synchronized cuts on every paired binary
word of length at most eight, and the finite cyclic-factor invariant with exact
mass `1/p`.

These finite checks support the stated helper claims. They do not prove
bi-infinite extendability, uniqueness of a maximal-entropy measure, or Bernoulli
measure isomorphism. Those conclusions require general arguments and external
ergodic-theory theorems, including the weak-Bernoulli isomorphism criterion used
by Wu and Zhang (arXiv:2611.03758v1). No Lean settlement is supplied here.
