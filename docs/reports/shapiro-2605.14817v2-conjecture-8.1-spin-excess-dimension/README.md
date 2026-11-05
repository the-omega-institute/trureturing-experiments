# Spin factors exceed the proposed component codimension

Source: Boris Shapiro, “Reducibility of spectral curves of finite Jacobi pencils”, arXiv:2605.14817v2 (13 August 2026), Sections 1, 3, 7.1 and Conjecture 8.1. Corollary 6.6 is the source of the spin factorization; its identity is prior work.

Run with Python 3 and SymPy:

```sh
python3 check.py
```

`output.txt` contains the complete output. All algebra uses exact integer polynomial coefficients (with symbolic indeterminates for the affine deformation). The independent determinant check uses a tridiagonal matrix with upper entries `w*c_i` and lower entries `w`, whose adjacent-entry products equal those of the symmetric source matrix with `b_i^2=c_i`.

## Result and scope

The check reproduces Corollary 6.6 for every integer size from 2 through 24, including the requested symbolic cases 4 and 6. It separately checks the determinant in sizes 4 and 6, and the proposed grouped-factor membership for sizes 16, 20, 24, 28 and 32. For size 16, all 16 diagonal entries are distinct and all 15 couplings are positive; the factor selecting source indices 5 through 12 has degree 8.

The mathematical obstruction is a dimension bound, not a numerical dimension estimate. The nonempty subset stratum lies in the 31-dimensional principal open subset of complex affine 31-space. Every nonempty irreducible component has nonnegative Krull dimension, so its codimension is at most 31. Conjecture 8.1 requires codimension 32.

For every integer m at least 4, use size n=4m, diagonal a_i=4m+1-2i, coupling c_i=i(4m-i), and subset S={m+1,...,3m}. Group the source's first m quadratic factors. Their specialization is exactly the product of (lambda+a_i) over S. Thus the stratum is nonempty, while its proposed codimension 2m² exceeds the ambient dimension 8m−1 by 2m(m−4)+1. The program verifies that polynomial identity symbolically; the nonnegativity for m≥4 is the elementary integer inequality displayed here. The general factorization is supplied by the source's Corollary 6.6, rather than inferred from the finite checks.

A stronger geometric bound comes from the affine spin family a_i=u+v(n+1−2i), c_i=zi(n−i), where v and z are nonzero complex numbers. Its quadratic factors are (lambda+u)²−(2r−1)²(v²+zt). Its inverse parameters are u=(a_1+a_n)/2, v=(a_1−a_n)/(2(n−1)), z=c_1/(n−1). This is a three-dimensional irreducible family in the same stratum, so some component at n=16 has codimension at most 28. The program checks the symbolic identity and inverse in sizes 4 and 6; the all-size deformation follows algebraically from the source identity by scaling and polynomial descent. This stronger bound is optional for the proposed settlement.

The computation checks polynomial equalities, distinctness, positivity, and integer dimension arithmetic. It does not compute irreducible components or prove the geometric statements in Lean. The preregistration proposes refuting only Conjecture 8.1; it does not settle the paper's reversal, degree-eight, or minimum-codimension conjectures.
