# Lang's A130187 odd-part numerator conjecture

Run `python3 check.py` in this directory. The checker needs Python 3 and its standard
library, uses integers and `fractions.Fraction`, and ends with `ALL_CHECKS_OK` on success.
`check.out` is its complete output.

The formal identity `(exp(g))' = g' exp(g)`, with `g = X^2/4`, computes the EGF
coefficients independently of the candidate double-factorial product. All 30 displayed
`r(n)` values and all 40 displayed `z(n)` values in Lang's linked note agree. The note's
label says `z(n), n=0..40`, while the displayed list has 40 entries, through index 39.
The EGF also gives `z(40)=0`.

All 1001 indices `0 <= n <= 1000` satisfy the proposed numerator formula, its denominator
formula, the coefficient bridge, and the disclosed LODA gcd formula. The checker also
compares every one of the 404 published b-file rows. These are finite exact arithmetic
checks; they do not establish the all-index assertion or constitute Lean verification.

Sources and attribution:

- The Online Encyclopedia of Integer Sequences, [A130187](https://oeis.org/A130187),
  revision 19; conjecture by Wolfdieter Lang, Jan 04 2013.
- Wolfdieter Lang, [Rationals and z-sequence](https://oeis.org/A130187/a130187.pdf),
  Jun 01 2007: the `SOURCE_R` and `SOURCE_Z` numerical fixtures.
- G. C. Greubel, [A130187 b-file](https://oeis.org/A130187/b130187.txt):
  `b130187.txt` is the unmodified downloaded snapshot; its SHA256 is printed by the checker.
- [LODA A130187](https://github.com/loda-lang/loda-programs/blob/main/oeis/130/A130187.asm)
  already implements the double-factorial/gcd generator. Formula discovery is not claimed.

The OEIS data and quoted fixtures are distributed under
[CC-BY-SA-4.0](https://creativecommons.org/licenses/by-sa/4.0/), as specified by the
[OEIS End-User License Agreement](https://oeis.org/LICENSE).
