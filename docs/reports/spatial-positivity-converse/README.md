# Spatial positivity and modular coefficients

The proposed counterexample concerns Buisine–Casini–Martinek, *Finding CFTs from spatial positivity of partition functions*, arXiv:2610.30234v1, the assumed equivalence after (2.16), repeated after (5.2).

Run `python3 check.py`. Only Python's standard library is needed. All evidence uses integers, rational numbers, or exact pairs representing Q(sqrt(2)); no floating point is used. The output is in [check.out](check.out).

The fidelity check reproduces the Ising projector decompositions in §3.2.1 (TeX labels eq:TTDLIsing and eq:projectors_verlindes_ising) and precisely the three nonnegative combinations of simple TDLs stated there.

For the candidate, the program checks every entry of the 8-by-8 rational transformed matrix, character equality under exchanging distinct labels 101 and 110, and the unique vacuum. For r=0,…,6 it checks the full transformed matrix against C tensor I and the exact partition-function weight polynomial x(x+y)^(r+2). The appended-bit Fourier identity records the algebraic recurrence used for arbitrary r. The first 25 character-series coefficients are checked exactly after extracting their fractional leading powers.

These are algebraic and finite-series checks. Absolute convergence, theta/eta modular transformations, the realization as genuine SU(2)_1 characters, and spatial positivity at every finite complex upper-half-plane test list require the analytic argument in the preregistration and a subsequent formal proof. The computation does not establish those analytic statements.

Source formulas and modular data are attributed to the paper and to Munizzi–Schnitzer, arXiv:2510.15067v3, §§2.4–2.5. The check is independently written; no source code or source-paper files are redistributed.
