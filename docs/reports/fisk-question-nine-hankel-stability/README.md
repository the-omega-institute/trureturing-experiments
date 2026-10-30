# Fisk Question 9: coefficient-Hankel stability

Source: Steve Fisk, *Questions about determinants and polynomials*, arXiv:0808.1850v1, Question 9 and matrix (4), Definition 1 and Remarks 1 and 8. Source statements are cited, not redistributed.

Run with Python 3 and SymPy:

    python3 docs/reports/fisk-question-nine-hankel-stability/check.py

The program uses exact symbolic algebra, integers and rational numbers. `output.txt` is its complete output; success ends with `ALL_CHECKS_OK`. No approximate root computation is used.

Fidelity checks reproduce Remark 1's degree-two coefficient transform and discriminant, and the d=1 assertion of Remark 8 for (x+y+1)^4.

For every real a>0, the candidate input is f_a(x,y)=(ax+y+1)^2(x+y+a)^2. The program checks positivity of its 15 supported coefficients and the positive imaginary parts of its affine factors. It computes the prescribed i=0,d=3 determinant by an independent Leibniz expansion and SymPy's determinant, obtaining

    H_a = 5A^4 + 8A^3B + 9A^2B^2 + 8AB^3 + 5B^4,
    A=ax+1, B=x+a.

Exact quadratic factorizations, discriminants, resultants and root transport identities give the stability condition a+1/a <= 16+4 sqrt(21). At a=1 the determinant is 35(x+1)^4; otherwise it has four distinct roots on the unit circle. Every a>=35, and every 0<a<=1/35, gives exactly two roots with positive real part. Fifty integer determinant recomputations and five rational sign checks supplement the symbolic checks.

The bottom-right entry of source matrix (4) is printed as f_(2d), without i. The witness has i=0, so it satisfies both the printed matrix and the fully shifted Hankel interpretation.

These are exact algebra checks supporting a preregistration candidate. They are not a Lean proof or an exhaustive literature clearance.
