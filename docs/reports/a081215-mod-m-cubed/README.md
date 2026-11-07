# Englander A081215 Conjecture XX

`check.py` uses Python 3 standard-library integers and `fractions.Fraction` only.
Run `python3 check.py`; `output.txt` is its complete output. Success requires
exit 0 and the final line `ALL_CHECKS_OK`.

Source: Mathew Englander, *Comments on OEIS A081215*, October 19, 2020,
<https://oeis.org/A081215/a081215.pdf>, definition on physical page 1 and
Conjecture XX on physical page 69 (zero-based index 68). The 72-page PDF's
SHA-256 is `74b6cccddfbb6b9ba54e501eed0959f0324d1f60c7e661dfc4b344c1c7ba0f43`.
The sequence is `a(n)=(n^(n+1)+(-1)^n)/(n+1)^2`, indexed from zero.
The conjecture asks whether `a(2m+1) mod m^3 = m` for every integer `m>=2`.

The fidelity checks reproduce both page-69 examples, the 19 displayed OEIS
values, and the proved source Theorem I for `n=1..2001`. The settlement checks
use both exact quotients and an independent modular-power calculation for
all `m=2..1000`. The binomial-tail and denominator-cancellation identities are
checked for all `m=2..128`; the inverse polynomial identity is also checked by
exact coefficient multiplication. Integer conversion occurs only after the
rational denominator is checked to equal one. No inverse of the even factor
four is taken modulo `m^3`.

These finite checks are evidence for a proposed universal proof, not a proof
of the universal assertion. No Lean proof is included.
