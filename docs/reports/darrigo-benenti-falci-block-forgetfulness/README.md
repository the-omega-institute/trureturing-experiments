# D’Arrigo–Benenti–Falci block-forgetfulness check

`check.py` uses SymPy exact expressions and rational arithmetic.  It verifies the one-use Gaussian dephasing fidelity obtained from the source’s factor `g`, the Lorentzian specialization of the covariance tail, and the covariance-cut boundary sum for every `1 ≤ N ≤ 4`, `2 ≤ M ≤ 5`, and `0 ≤ L ≤ 4`.  The exact choice `τ = log(3)`, `τ_p = τ_c = 1` gives `r = 1/3` and `τ > τ_p`, so the separated-window formula is valid.  It also checks the displayed `N=1,M=2` Lorentzian bound factor from Eq. (38).

No floating point values are used as evidence.  The generated output is stored in `check.out`.
