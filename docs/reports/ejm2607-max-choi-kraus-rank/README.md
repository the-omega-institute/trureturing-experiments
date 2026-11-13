# Maximum Choi–Kraus rank with a unit

Exact rational checks accompanying a proposed answer to Evert–Jury–McCullough,
*Completely Positive Matrix Products*, arXiv:2607.13251v1, Section 5.3, Question 1.
Source: https://arxiv.org/abs/2607.13251v1.

Run `python3 docs/reports/ejm2607-max-choi-kraus-rank/check.py` from the repository
root. Python 3.12 or newer suffices; the program uses only the standard library.
The captured output is in `output.txt` and ends with `ALL_CHECKS_OK`.

The fidelity check reproduces Example 5.10: its doubled-diagonal weighted Schur
product, both unit identities, Choi formula, and an indefinite unit. The proposed
settlement is checked for n = 1, 2, 3. For n ≥ 2 the proposed product is

Φ(A,B) = (I+J) ∘ A ∘ B + Tr(HA)Tr(HB)I,
H = (2n−1)I−J, T = J−I/2.

Here J is the all-ones matrix and ∘ denotes entrywise multiplication. The script
evaluates every Choi entry independently on matrix units, verifies both unit
identities on a basis, and checks the Gram identity for K = nI−J, the Choi
Kronecker expansion, and the bound C ≥ (n−1)²I by exact symmetric elimination.
Independent row elimination cross-checks ranks and determinants. All arithmetic
uses integers and `fractions.Fraction`; no floating point is used as evidence.

These finite checks support the all-dimension paper argument. They prove neither
the infinite family nor a Lean theorem. The source's Question 1 and the full
preregistration statement include every positive integer n, with output size
m = n, and impose no associativity, causality, trace condition, or positivity
condition on the unit.
