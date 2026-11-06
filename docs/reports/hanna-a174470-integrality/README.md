# Hanna A174470: exact coefficient and identity checks

The target is Paul D. Hanna's OEIS A174470 conjecture that every coefficient of

\[
\exp\left(\sum_{n\ge1} b(n)3^n x^n/n\right)
\]

is an integer. Here b(0)=0 and b(n)=1 exactly when the 3-adic valuation of n is even.
This is the Thue/Roth sequence A014578, distinct from the Thue–Morse sequence.

Run with Python 3, using only the standard library:

    python3 docs/reports/hanna-a174470-integrality/check.py

The stored result is [output.txt](output.txt). All arithmetic used as evidence is
integer or `fractions.Fraction` arithmetic. The checks cover:

- Fidelity to the A014578 substitution morphism at indices 1 through 2187, and its
  recurrence and floor-sum formula for n=1 through 10000.
- All 27 A174470 coefficients displayed in its source entry, at indices 0 through 26.
- Integer denominators equal to 1 through degree 120, computed from the formal
  exponential's differential recurrence over the rationals.
- Agreement through degree 120 with the locally finite fractional Euler product
  whose factor at d=3^k is (1-3^d*x^d)^((-1)^(k+1)/d), including its explicit
  factorial coefficient formula.
- The proposed Dold congruences for every prime p at most 101 and m=1 through 300.
- The sharp-extension prediction for integer bases B=-12 through 12, through
  degree 45, including the obstruction coefficient 2*B^3/3.

These are finite checks and do not establish the all-index conjecture. The
planned settlement uses a general prime-power divisibility argument and the
existing Dold exponential-integrality theorem; finite enumeration is only a
check of definitions and identities.

Sources: [A174470](https://oeis.org/A174470/internal), revision #2,
last updated March 30, 2012; originally submitted March 28, 2010;
[A014578](https://oeis.org/A014578/internal), revision #68,
last updated October 14, 2026. Both entries were retrieved October 28, 2026.
The displayed coefficient list is transcribed into the program for reproducible
fidelity checks. OEIS data are licensed under
[CC BY-SA 4.0](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement).
