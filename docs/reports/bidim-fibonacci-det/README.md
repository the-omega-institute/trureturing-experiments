# Bidimensional Fibonacci determinant checks

Source: Eudes Antonio Costa and Paula M.M.C. Catarino, *On Matrices with
Bidimensional Fibonacci Numbers*, Annales Mathematicae Silesianae 40(1),
59-76, DOI [10.2478/amsil-2025-0014](https://doi.org/10.2478/amsil-2025-0014).
The target is the unnumbered Conjecture on p. 73.

Run `python3 docs/reports/bidim-fibonacci-det/check.py` from the repository
root. Python 3's standard library suffices. `out.txt` is the resulting output.
Gaussian integers are pairs of Python integers; no floating point is used.

The program constructs the sequence from Definition 2.1, compares it with
Lemma 2.5 on 10,609 entries, and reproduces Example 4.1 and the formula in
the proof of Theorem 4.3. For every pair `0 <= m,n < 100`, it checks the
matrix determinant and the identity

`Im(det F_(m,n)) = (-1)^m T_n + (-1)^n T_m`,

where `T_r = fib(r+1)*fib(r+2)`. All 10,000 imaginary parts are nonzero.
The finite checks support a proof outline; they do not prove the universal
conjecture, and this report contains no Lean proof.

The program also records a printed discrepancy in Proposition 4.6:
at `(m,n)=(1,0)` its displayed final formula gives `-7-2i`, while
Definition 3.1 and Example 4.1 give `-7+i`. The checks use Definition 3.1.
