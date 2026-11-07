# OEIS A126364 adjacent-digit boundary check

`check.py` uses only exact integer dynamic programming, exact `fractions.Fraction`, and SymPy's exact series algebra. It reproduces the first 27 terms of OEIS A126364 and A025565, then verifies

\[
a_n-3a_{n-1}+b_n=2(U_{n-2}-V_{n-2}),
\]

for every `2 <= n <= 200`, where `U_r` counts nonnegative Motzkin walks and `V_r` counts those additionally bounded by 9. The same run checks equality for `2 <= n <= 11`, strict failure for `n >= 12`, and the exact `n=12` defect of 2. The all-`n` conclusion in the preregistration is supplied by the reflection and endpoint-extension bijections; this script is the finite exact regression check.
