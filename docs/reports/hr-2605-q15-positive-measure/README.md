# Positive measure for constant-section parameters

Source: J. Haddad and D. Ryabogin, *On convex bodies with constant non-central sections*, arXiv:2605.00299v1, Question 15 and Definition 8.

Run `python3 check.py`. The program uses integers and `fractions.Fraction` throughout, with outward rounding to a rational grid. It encloses π by Machin's formula, square roots by integer square roots, and sine, cosine and arctangent by Taylor series with explicit remainder bounds. The literal source sine quotient is used for F₊, with an independent overlap check against its cosine form. There is no floating point evidence or external Python dependency.

The checks reproduce the source's continued fraction and convergent denominators for γ₀=(4+√2)/14 and its inequalities C(A₀,k) for k=2,…,6. For γ*=(40+√2)/94, they certify R*>5, the proposed upper and lower bounds at k=2,3,4, and all inequalities through k=150. They also test the two proposed Diophantine bounds in all 2,346 eligible small cases described in the output, the exponential majorant starting at J=160, and finite instances of the cylinder tail estimate.

The expected affirmative settlement uses a general argument: two odd denominators among adjacent convergents and their difference control the second smallest odd-multiple error; best approximation bounds the next error from below; an exact seed and continuity handle a finite prefix; summable digit exclusions leave a positive-measure subset of a continued-fraction cylinder. The checks certify the finite seed cases and test intermediate claims. They do not prove the infinite assertion, and no Lean proof is included.

The stored output is `check-output.txt`; a successful run ends with `ALL_CHECKS_OK`.
