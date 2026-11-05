# Seven-data positive Girard–Waring coefficient checks

Run `python3 docs/reports/wojnar-a301424-seven-data-gf/check.py` from the repository root. Python 3 with its standard library is sufficient. The recorded output is in `output.txt`; successful execution ends with `ALL_CHECKS_OK`.

The coefficient formula is equation (2) of Wojnar, Wojnar and Brin, [arXiv:1706.08381v2](https://arxiv.org/abs/1706.08381v2). The definitions and 21 seven-data terms come from [OEIS A301417](https://oeis.org/A301417) revision 42 and [A301424](https://oeis.org/A301424) revision 26. Each multiplicity vector represents one monomial in independent normalized elementary symmetric variables. Positivity is tested before evaluating those variables on a data set.

The program reproduces the paper's three-data quadratic and cubic polynomials and all seven six-data table masses. For seven data, at every degree 1 through 27 it compares the complete coefficient dictionary from the factorial formula with an independently implemented Newton recurrence, and compares strict positive mass with exact formal division of the conjectured rational function. All 21 published terms agree. Signed mass is 1; the magnitude of the negative mass is one less than the positive mass. Four exact data families also verify evaluation against the actual power means, including repeated, zero, negative and rational data.

At degree 27 there are 1175 monomials. Positive mass is 24203797992661339086776056; negative magnitude is 24203797992661339086776055. The expanded generating-function denominator, in increasing degree order, is `[1, -8, -14, -14, 0, 14, 14, 6, 1]`.

Only integers and exact rational arithmetic are used. These are finite source-fidelity checks and tests of the proposed settlement. They establish no identity at untested degrees and contain no Lean proof.
