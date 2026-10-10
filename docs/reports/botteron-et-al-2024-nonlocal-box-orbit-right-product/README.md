# Botteron–Broadbent–Chhaibi–Nechita–Pellegrini: the right product maximizes CHSH over the whole orbit

`check.py` supports a proof of Conjecture 15 ("Dyck paths") of arXiv:2312.00725v3 (Quantum 8, 1402 (2024)): for the Brunner–Skrzypczyk wiring W_BS on conv{PR, SR, I}, the right product P^{⊠k} attains the highest CHSH value over all parenthesizations of k copies of P.

- From the source's multiplication table, in coordinates A(u,a) = a PR + (u−a) SR + (1−u) I: (u,a)⋆(v,b) = (uv, av + (1+u−2a)b/4) and CHSH = 1/2 + (u+a)/4 (sympy).
- The associator identity a((A⋆B)⋆C) − a(A⋆(B⋆C)) = c(2(1−v)(2−t) + (1+v)t)/16 with t = 1+u−2a (sympy).
- All parenthesizations for k ≤ 9 (200 random boxes) and k ≤ 6 (2000 random boxes) of conv{PR,SR,I}, exact rationals: no parenthesization exceeds the right product.
- The multiplication table follows from the source's box-product formula with the 32-entry vector of W_BS (Section 4).

Run: `python3 check.py` with sympy (exit 0; lines `product and CHSH formulas: True`, `associator identity: True`, `parenthesization violations: 0`, `multiplication table from W_BS: True`, `bad= 0`).
