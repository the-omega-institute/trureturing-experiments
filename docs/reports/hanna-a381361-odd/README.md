# OEIS A381361: exact checks of the all-odd observation

Paul D. Hanna's [OEIS A381361](https://oeis.org/A381361), published revision 25 (April 23, 2026), defines the exponential generating series by

    A(x) = exp(Integral abs(1/A(x)^2) dx).

The first COMMENTS line says “All terms appear to be odd.” This report concerns that line. The separate modulo-five conjecture is outside this report's target.

Run with Python 3 and the standard library:

    python3 check.py

The integer binomial-convolution recurrence reproduces every source term at indices 0 through 500. All 501 terms are odd, and the positive-degree factorial-normalized coefficients of A squared and its inverse are even through degree 500. A separate implementation with ordinary rational coefficients evaluates the literal map using ordinary multiplication, unit inversion, a zero-constant primitive, and the finite Taylor exponential. Its fixed point through degree 18 agrees with the integer recurrence. Every comparison is exact; these are finite checks, not a proof of the infinite assertion.

`output.txt` is the program's complete output, ending with `ALL_CHECKS_OK`. `source_terms.tsv` is the unchanged [source b-file](https://oeis.org/A381361/b381361.txt), attributed to Paul D. Hanna and the OEIS contributors; SHA-256: `93779794cd57e856400e4156f8eb50cc80971ebed6549b7c2f6647ce6f71dfd4`. The source data retain [CC BY-SA 4.0](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement). The independently authored check and this report are also available under CC BY-SA 4.0.
