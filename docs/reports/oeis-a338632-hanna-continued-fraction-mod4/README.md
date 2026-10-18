# OEIS A338632 (Hanna): a(n) ≡ 2 (mod 4) for n > 1

Entry: "G.f. A(x) satisfies: 1 = A(x) - x/(A(x) - 3*x/(A(x) - 5*x/(A(x) - 7*x/(A(x) - ...)))), a continued fraction relation." Conjecture: "a(n) = 2 (mod 4) for n > 1".

Route: reduce modulo 4. Since 2(j + 2) + 1 ≡ 2j + 1 (mod 4), the tails of the continued fraction are 2-periodic.
- **The quadratic.** With U = T₀ = 1 and V = T₁ = A − 3x, the relation becomes (A − 1)(A − 3x) = x.
- **The closed form.** B = (1 + x²)/(1 − x) satisfies the same quadratic modulo 4.
- **Uniqueness.** The factor A + B − 1 − 3x is a unit, so A ≡ B = 1 + x + 2Σ_{n≥2} xⁿ (mod 4).

## Program

`check.py` (standard library only). Run `python3 check.py`; it prints `RESULT PASS` when every check holds. It computes A to order 40 exactly as the fixed point of A ↦ 1 + x/T₁(A), and checks that:
- every coefficient is an integer and agrees with the OEIS data;
- the defining relation holds to order 40;
- a(n) ≡ 2 (mod 4) for 2 ≤ n ≤ 40;
- A ≡ (1 + x²)/(1 − x) (mod 4).

It also reports that A² ≡ A + x (mod 3), which bears on the entry's separate mod-3 conjectures.

Reading: `run.log` (exit 0, `RESULT PASS`).
