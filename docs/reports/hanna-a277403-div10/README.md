# Hanna A277403: exact divisibility checks

Source: [Paul D. Hanna, OEIS A277403](https://oeis.org/A277403), revision #21, offset 1. The COMMENTS conjecture asks whether `a(n)` is divisible by 10 for every `n > 2`. Source terms are attributed to the OEIS under its [End-User License Agreement](https://oeis.org/LICENSE).

Run with Python 3 and its standard library:

```sh
python3 docs/reports/hanna-a277403-div10/check.py
```

`output.txt` records a successful run. Integer Bell recursion reproduces all 23 displayed source terms and finds no divisibility failures through index 60. Independent ordinary power series substitution with `fractions.Fraction` agrees through index 24. Recursions modulo 2, 5 and 10 agree with the integer reduction through index 60.

For `P=x+x^2` and `J=x^2/2+x^3/3`, the complete residual `P(x-J)-x-J` has exponential coefficients `-10,-10,40,80` in degrees 3, 4, 5, 6 and zero elsewhere. This identity is checked over the rationals before reducing its integral exponential coefficients. No denominator is inverted modulo 10.

These are finite corroboration and an exact polynomial identity. The all-index integral construction and uniqueness modulo 10 remain proof obligations; this report contains no Lean proof of the conjecture.
