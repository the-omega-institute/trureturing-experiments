# Lin, Conjecture 4.3: the additive majorization fails at p = 4

Check for the open-problem preregistration refuting Conjecture 4.3 of M. Lin, "On a determinantal inequality arising from diffusion tensor imaging", Communications in Contemporary Mathematics 19 (2017), 1650044 (arXiv:1604.04141).

Run: `python3 check.py` (Python 3 with sympy and numpy; the certificate uses exact rationals).

With A = diag(3, 1, 1/2) and B = [[1/2, 1/4, 3/4], [1/4, 3/4, 0], [3/4, 0, 5/2]], put L = A² + (AB²A)² = A² + |BA|⁴ and R = A² + (BA²B)² = A² + |AB|⁴. The program checks, exactly:
- A and B are positive definite (leading principal minors);
- tr L = tr R;
- L − (7/5)I is positive definite, so the smallest eigenvalue of L exceeds 7/5;
- det(R − (7/5)I) < 0, so the smallest eigenvalue of R is below 7/5.

With equal traces, the sum of the two largest eigenvalues of R therefore exceeds that of L, so λ(L) does not majorize λ(R). In floating point it also prints both spectra, checks that the same pair satisfies Lin–Wolkowicz (2012), Corollary 2.3, at k = 2 (which compares A² + (AB)²(BA)² with A² + (BA)²(AB)², a different pair of matrices), and lists the exponents p on a grid of step 1/4 in [1/4, 12] at which this pair violates the majorization.

`run.log` holds the output of the run cited by the preregistration.
