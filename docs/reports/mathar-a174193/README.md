# Mathar's OEIS A174193 recurrence

`check.py` checks the order-four polynomial recurrence conjectured by R. J. Mathar in the FORMULA field of [OEIS A174193](https://oeis.org/A174193/internal). The object is the expansion of the displayed algebraic generating function, with the square-root branch of constant coefficient one. It does not assume or settle the identification of that generating function with the permutation-counting series.

Run `python3 docs/reports/mathar-a174193/check.py` from the repository root. Dependencies: Python 3 and SymPy. `output.txt` contains the successful output and dependency versions. Every calculation uses integers, fractions, or exact polynomial algebra; no floating point is used as evidence.

The checks reproduce all 30 coefficients displayed by OEIS, compare the hook-length sum in Section 2 of [Almost avoiding permutations, arXiv:2007.15111v1](https://arxiv.org/abs/2007.15111v1) for indices zero through 25, and independently enumerate permutations of lengths zero through seven under the paper's almost-321-avoidance definition. They check the proposed recurrence at indices four through 200. These finite checks are corroboration.

The symbolic part verifies the shifted Euler-operator conversion and two explicit polynomial certificates used in the universal differential-identity argument. It does not provide a Lean proof or claim a kernel-checked settlement. The source TeX dates the background article June 22, 2009; arXiv v1 was submitted July 29, 2020. The recurrence itself is OEIS-only, attributed June 14, 2016, and is absent from the background TeX.

Source attribution: R. J. Mathar and the OEIS Foundation, A174193 (OEIS content under [CC BY-SA 4.0](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement)); Robert Brignall, Shalosh B. Ekhad, Rebecca Smith and Vincent Vatter, *Discrete Mathematics* 309 (2009), 6626–6631, [DOI](https://doi.org/10.1016/j.disc.2009.06.027). This report contains independently written checking code and generated results.
