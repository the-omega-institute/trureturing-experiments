# Pentarnion nuclei: exact independent checks

Source: Alireza Poshtkohi and Susanne Pumplün, *From Quaternions to
Pentarnions I: Algebra, geometry, and a spin model*, arXiv:2610.31485v1
(22 October 2026), Sections 1, 2 and 7, Table 2 and Theorem 3.9.
<https://arxiv.org/abs/2610.31485v1>

Run with Python 3 and Sympy 1.14.0:

```sh
python3 check.py > output.txt
```

The program uses exact integers, rationals, symbolic polynomial identities
and exact real algebraic radicals. It uses no floating point. Its
multiplication tensor independently transcribes the source basis products.

It reproduces Theorem 3.9's formula for the normalized associator
[x, ell, ell] and its coefficient determinant 54. It checks scalar
inclusion, conjugation reversing products, the left and middle nucleus
elimination identities, and the polynomial identities used to show that
the three mixed coefficients cannot vanish.

For every tuple with a, b, c, d in {-2, -1, 1, 2} and gamma in
{-2, -1, 0, 1, 2}, it applies the source's admissibility conditions exactly.
There are 126 admissible tuples, including 46 with indefinite diagonal
norm and 22 with gamma = 0; the remaining 1154 tuples are outside the
source domain. For each admissible tuple, a nonzero four by four minor of
actual associator equations forces the pure part of a left or middle
nuclear element to vanish. Conjugation gives the right nucleus. Unknown
answers from Sympy's exact zero predicate fail the check.

The output ends with ALL_CHECKS_OK. The finite grid is evidence supporting
the proposed general argument, not a proof for all real parameters. The
automorphism-group conjunct of the source's question is outside this check.
No Lean proof of the claim or its negation is included.

The source is cited, not redistributed. The program is original code
written for this experiment, under this repository's license.
