# Bala's fractional factorial ratios (OEIS A347854–A347858)

Checks for the open-problem preregistrations on the integrality of the five sequences of P. Bala in OEIS A347854–A347858, and of the auxiliary b(n), c(n) in A347855. Fractional factorials are x! := Γ(1 + x).

Run: `python3 check.py 300` (Python 3, standard library only; exact rational arithmetic).

The program removes the Gamma factors exactly:
- Γ(1 + t/2)/√π = t!/(2^t ((t − 1)/2)!) for odd t, and √π cancels in every ratio;
- Γ(1 + 4n/3)/Γ(1 + n/3) = ∏_{j=1..n} (n/3 + j).

It checks:
- the published data;
- integrality for n ≤ 300;
- a(3m+1) = 3^(3m+1)·b(m) and a(3m+2) = 3^(3m+2)·c(m);
- nonnegativity of the ten step functions H, K and T_r of the proofs on the grid j/360 over one period. Every breakpoint lies on this grid, and the functions are right-continuous and constant between breakpoints, so this check proves nonnegativity everywhere.

`run.log` holds the output of the run cited by the preregistrations.
