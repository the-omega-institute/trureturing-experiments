# Hanna A231899: congruences modulo 5 and 100

Run `python3 check.py` with Python 3; no third-party packages are required. The stored output ends with `ALL_CHECKS_OK`.

The source is [Paul D. Hanna, OEIS A231899](https://oeis.org/A231899), approved revision 19, 2013-12-03. Its defining equation is A(0)=1 and A'(x)=A(x(A'(x))^3). The nominated conjectures say that a(n)=n![x^n]A(x) is divisible by 5 for n≥5 and by 100 for n≥10. The embedded source data are its 19 displayed terms and b-file entries 0 through 40, from [the b-file](https://oeis.org/A231899/b231899.txt). OEIS content remains subject to the [OEIS End-User License Agreement](https://oeis.org/LICENSE); new code is MIT licensed.

The program constructs the sequence independently by the integer partial-Bell-polynomial recurrence and by the source's Lagrange coefficient formula. They agree at indices 0 through 40 and reproduce every embedded source value. Both nominated congruences hold throughout this computed range.

For F_d=Σ(j=0..d) a(j)x^j/j! it expands R_d=F_d'−F_d(x(F_d')^3) using exact rational arithmetic. It checks all 41 coefficients for d=4 and all 226 coefficients for d=9, after factorial normalization, for integrality and divisibility by 5 and 100 respectively. It also checks the displayed low residual coefficients, denominator conditions, the weighted 2-adic coefficient bounds of F_9, F_9' and R_9, and the factorial bound at indices 16 through 10000. Zero coefficients are handled separately from the integer-valued valuation function.

The program verifies three obstructions to direct application of previously published general congruence theorems: [y^4]A(log(1+y))=23/8, w_4=27/2 and w_17=114101705219927291073/2242240 for w_m=m[x^m]log A. The last denominator is divisible by 5.

These are finite supplemental checks. An unbounded settlement requires the triangular transfer argument and the all-index denominator and valuation estimates described in the preregistration; neither this output nor this script is a Lean proof or a completed settlement.
