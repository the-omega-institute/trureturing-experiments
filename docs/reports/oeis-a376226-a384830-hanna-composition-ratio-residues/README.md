# OEIS A376226 and A384830 (Hanna): residues modulo 2 and 3

Entries:
- **A376226.** U(x) = U(x³ + 6xU³)/U(x² + 4xU²), with conjectures (C1) a(n) ≡ 1 (mod 3) iff n = 3^k and (C2) a(n) ≡ 2 (mod 3) iff 2n = 3^i + 3^j with i < j.
- **A384830.** U(x) = U(x³ − 3xU³)/U(x² − 2xU²), with conjectures (C.1) a(n) is odd iff n = 2^k, and (C.2), (C.3) as above.

## Route (GPT Pro r386q; every step checked by the orchestrator)

- **The recursion.** It has multiplier −1, so the solution is integral and unique over ℤ, F₂ and F₃.
- **Modulo 2.** The A384830 equation is solved by S = Σx^{2^k}.
- **Modulo 3.** Both equations become B(x³) = B·B(x² + xB²). This is solved by the inverse B of y + y² + y³, and B(t²) = (Σt^{3^k})² gives the residues.

## Program

`check.py` (standard library only; reads the saved entry texts). Run `python3 check.py`; it prints `RESULT PASS` when every check below holds.
- **The recursion.** Its multiplier is −1 at every step.
- **The entries.** The first 12 terms of both entries match.
- **The conjectures.** They hold for n ≤ 70.
- **The reductions.** Both sequences agree modulo 3 with the inverse of y + y² + y³, and A384830 agrees modulo 2 with Σx^{2^k}.

Reading: `run.log` (exit 0, `RESULT PASS`).
