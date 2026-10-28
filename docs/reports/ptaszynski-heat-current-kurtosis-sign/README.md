# Ptaszyński thermal heat-current kurtosis sign

`check.py` is an exact-arithmetic reproduction and counterexample check for Eq. (19) of K. Ptaszyński, arXiv:2204.08421v4, *Bounds on skewness and kurtosis of steady-state currents*. It builds the heat-tilted generators and determinants with SymPy, solves the implicit dominant-eigenvalue series through fourth order, reproduces a positive two-state thermal case, and verifies the six-state locally detailed-balanced family with `c₂ > 0` and `c₄ < 0` at `R = 8` (and the supplied polynomial sign certificate for every integer `R ≥ 8`). A second witness has six distinct energies.

The program uses integers, SymPy exact algebra, and rational arithmetic only. It prints `ALL_CHECKS_OK` on success and exits with status 0.
