# Hanna A227051 modulo 10

Source: Paul D. Hanna, [OEIS A227051](https://oeis.org/A227051), revision 10.
The first FORMULA conjectures `a(n) == 3 (mod 10) for n>0`, for the even,
factorial-normalized coefficients of `1/(cos(x)-sin(x)*sinh(x))`.

Run `python3 check.py` with Python 3; only the standard library is required.
`output.txt` contains the exact output. The program uses integers and
`fractions.Fraction` throughout and makes no network calls.

It independently inverts the literal rational Taylor series through degree
240 and an integer binomial recurrence through index 120. It reproduces all
16 published DATA values, checks vanishing odd coefficients, integrality,
the proposed denominator formulas, the even and alternating binomial sums,
and the inverse residues modulo 2 and 5. Index 0 has value 1 and is excluded.
It also checks `6*[x^6]log(A)=379/60`, which fails the modulus-integral
logarithmic-weight hypothesis of Kallat's shift-congruence theorem at 2 and 5.

These finite checks support transcription and a proposed proof route; they
do not establish the unbounded conjecture. No Lean proof is included.
