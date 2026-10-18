# OEIS A380066 (Hanna): a(n) ≡ 2 (mod 4)

Entry: "G.f. A(x) satisfies 1/2 = Sum_{n=-oo..+oo} x^n*A(x)^n * (A(x)^n + x)^(n-1) * (A(x) + x^n)^(n-1)." Conjecture: "for n > 0, a(n) == 2 (mod 4)."

This is the case c = 1 of the family settled for c ≥ 2 in the frozen `D5.S1.Recurrence.Bilateral.BilateralProductModFour` (A381362, A381363).

## Route (GPT Pro r386q; checked by the orchestrator)

- **The negative summands.** The k = −m summand equals (xA)^{m²}/((1 + xA^m)^{m+1}(1 + x^mA)^{m+1}). So the equation is 1/2 = 1/((1 + x)(1 + A)) + S(A), with S(A) integral and contracting.
- **Modulo 4.** At the fixed point (1 + x)(1 + A) = 2/(1 − 2S) ≡ 2 (mod 4), so A ≡ (1 − x)/(1 + x) (mod 4).

## Program

`check.py` (standard library only; reads the saved entry text). Run `python3 check.py`; it prints `RESULT PASS` when every check below holds.
- **The rewrite.** It reproduces the first 10 terms of A380066, and with c = 2 the frozen A381362.
- **To order 18.** A is integral, a(n) ≡ 2 (mod 4), and A ≡ (1 − x)/(1 + x) (mod 4) coefficientwise.

Reading: `run.log` (exit 0, `RESULT PASS`).
