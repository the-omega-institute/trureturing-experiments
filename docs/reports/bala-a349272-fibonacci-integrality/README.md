# Bala's even-step Fibonacci polynomial integrality

The source is Peter Bala's second COMMENT in [OEIS A349272](https://oeis.org/A349272/internal), internal revision #13. The preceding COMMENT restricts the step to a positive even integer; the entry's offset is zero. Fibonacci polynomials have initial values F_0 = 0 and F_1 = 1 and recurrence F_(j+2) = x F_(j+1) + F_j.

Run from this repository's root:

```sh
python3 docs/reports/bala-a349272-fibonacci-integrality/check.py
```

Python 3's standard library suffices. Coefficient tuples use integer arithmetic; rational evaluations use `fractions.Fraction`. No floating point is used. The output is stored in `output.txt`; success requires exit 0 and final line `ALL_CHECKS_OK`.

The fidelity check reproduces the ten displayed A349272 data terms at indices 0 through 9 and independently checks its first FORMULA. The identity checks cover 169 determinant-one recurrence addition identities, 13 odd partial sums, 72 even-step reductions, and all 36 pairs with m in {2,4,6,8,10,12} and n in {0,...,5}. For those 36 pairs, direct polynomial multiplication and independent monic long division agree on an integral quotient. The largest numerator degree is 781. A further 130 integer checks cover x=1, every even m from 2 through 20, and n from 0 through 12.

The proposed quotient uses c = F_m(x), t = F_(m+1)(x) + F_(m-1)(x), and U_0(t)=0, U_1(t)=1, U_(k+2)(t)=t U_(k+1)(t)-U_k(t). For n=0 it is 1. For n>=1 it is

H = c^(2n) (U_(n+1)(t)-U_n(t)) product_(1<=k<=2n, k!=n+1) U_k(t).

The program checks D H = P, where D is the source's denominator and P its numerator. It also checks that D is monic and nonzero, although D(0)=0 in every tested case. The quotient is understood in the rational-function field, with no assertion of pointwise division at x=0. The excluded odd step m=1,n=1 gives 1/2 at x=1.

These finite checks support a preregistration and do not prove the universal conjecture. The intended positive settlement requires a general algebraic argument. The program is original project code; its license is the repository's LICENSE. The OEIS numerical terms are cited source data, subject to the [OEIS End-User License Agreement](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement).
