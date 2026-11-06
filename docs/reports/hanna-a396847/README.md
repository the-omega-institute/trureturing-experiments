# Exact checks for Hanna's A396847 boundary conjecture

Source: Paul D. Hanna, [OEIS A396847](https://oeis.org/A396847/internal), revision 18, author line July 12, 2026; entry timestamp July 15, 2026. OEIS sequence data and the polynomial example are attributed to Hanna. No source paper is associated with the entry.

Run `python3 docs/reports/hanna-a396847/check.py` with Python 3. The program uses only integers and `fractions.Fraction` as evidence and requires no third-party packages. Its captured output is `check.out`; success ends with `ALL_CHECKS_OK` and exit 0.

The fidelity check reproduces all 36 coefficients printed in the source and the full F₂ polynomial. Further finite checks compare coefficient stabilization through degree 40, polynomial evaluation with the scalar recursion, and the lower trapping invariant on 63 short increasing exponent lists. Directed integer rounding at scale 10³⁰ bounds the prime-exponent orbits starting at 356/1000 and 357/1000. The exact geometric-tail inequality at 87/100 supplies the lower-end trapping certificate; the other orbit crosses 1.

These are finite evidence for a proposed general threshold argument. They do not prove the infinite coefficient/evaluation bridge or the boundary theorem, and they do not certify all digits of the OEIS decimal. The radius bracket is conditional on that bridge and the elementary continuity/tail argument stated in the preregistration. The planned main result is an unbounded mathematical argument, not a positive finite enumeration.
