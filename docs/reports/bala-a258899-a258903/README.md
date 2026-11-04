# Bala's totient-period conjectures for A258899 and A258903

`check.py` is a Python 3 standard-library program using integers and `fractions.Fraction` throughout. Run it with `python3 docs/reports/bala-a258899-a258903/check.py`. The captured output is `check.out`; its last line is `ALL_CHECKS_OK`.

Sources: [A258899 internal entry](https://oeis.org/A258899/internal), revision 17, and [A258903 internal entry](https://oeis.org/A258903/internal), revision 12, both published 9 August 2025 and retrieved 30 October 2026. Their Peter Bala conjecture comments are dated 3 August and 8 August 2025. The finite divisor-factorial definitions are the FORMULA fields; both offsets are zero with a(0)=1. A258903's displayed modulo-7 list starts at index 1.

The program reproduces A258899 at indices 0 through 19, A258903 at indices 0 through 16, and both source residue lists. An independent expansion of the locally finite exponential generating function agrees through degree 40. It checks the sources' displayed periods through index 240.

For the family c(u,v,0)=1 and c(u,v,n)=sum over de=n of u^d v^e n!/(d!e!), it checks the proper-factor quotient/product identity and divisibility, and the tail formula c(u,v,n)=u v^n+v u^n modulo m. The parameter range is 1<=u,v<=4, 1<=m<=60, 2m+4<=n<=240. The totient shift is checked whenever the shifted index is at most 240. Integer-power checks also include base zero. Every finite factorial quotient is first computed as a rational and checked integral.

Applicability checks distinguish these targets from existing general theorems. The inverse-Stirling ordinary coefficients include -2/3 for A258899 and 21/2 for A258903, excluding integral substitution. Both integer inverse-Stirling transforms fail Danesh's Theorem 3.1(ii) at p=5,r=8. The unique logarithmic cycle weights include 23/3 and -3269/4, excluding Kallat's integer-cycle theorem. References: [Danesh](https://doi.org/10.33774/coe-2026-4d6jr), [Kallat's integral-substitution theorem](https://arxiv.org/abs/2609.07238v1), and [Kallat's weighted-permutation theorem](https://arxiv.org/abs/2610.15849v1).

These are finite checks and prior-art applicability diagnostics. They do not prove either unbounded conjecture. The proposed settling route is a general symbolic argument, with these checks used only as independent evidence.
