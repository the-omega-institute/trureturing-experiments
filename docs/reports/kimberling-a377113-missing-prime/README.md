# Kimberling A377113: missing-prime certificate

Source: [OEIS A377113](https://oeis.org/A377113), revision 10. For the expansion of
`(3 + sqrt(2) + sqrt(3))^n` in the basis `1, sqrt(2), sqrt(3), sqrt(6)`, the source
conjectures that every prime divides the constant coefficient infinitely often.

Run `python3 check.py`. It uses Python 3 integers and standard-library functions
only. The stored output is `check.out`.

The check reproduces all 23 published sequence terms by two exact methods:
multinomial expansion with radical reduction, and four-coordinate multiplication.
It checks the proposed prime-167 certificate by that transition, an independent
fourth-order recurrence, and four split-field power evaluations. The full state
returns after 166 steps and every constant coordinate on the cycle is nonzero.
Determinism then supplies the all-index periodicity argument; this report contains
no Lean proof.

The source's FORMULA initial indices and EXAMPLE exponent are inconsistent with
its defining expression and DATA. The displayed small-prime gap patterns also
fail direct checking: the computed patterns for 2, 3, 5, 7 are respectively
`[1]`, `[2]`, `[8,1,4,3,8]`, and `[12,10,26]` through index 1000. These finite
comparisons are source errata, distinct from the proposed prime-167 refutation
of the general every-prime assertion. They are printed in the output.

The finite certificate is computational evidence for a preregistration, not a
claim of completed kernel verification or of exhaustive literature novelty.
