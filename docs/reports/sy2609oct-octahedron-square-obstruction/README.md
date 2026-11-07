# Octahedron diameter-plus-square obstruction

This checks the unnumbered scalar conjecture in section 3 of Jakob Steininger
and Sergey Yurkevich, *The Nieuwland number of the Octahedron*,
[arXiv:2609.10788v1](https://arxiv.org/abs/2609.10788v1).
For the convex hull O of the six signed coordinate vectors, let
Φ_c(A) = diam(A) + c s(A), where s is the smallest enclosing-square side,
allowing rotations and translations. The conjecture fixes c = 23/10 and
asserts sup Φ_c(MO) ≤ (3√2/4) inf Φ_c(MO) over MMᵀ = I₂.

Run `python3 check.py` with Python 3 and SymPy 1.14.0. The stored output is
`check.out`; its last line is `ALL_CHECKS_OK`. The program uses exact SymPy
algebraic expressions and outward rational intervals built with `Fraction`
and integer square roots. It uses no floating point or sampled angles.

The fidelity check reproduces the source's displayed M₁ vertices and the
side √2 of its M₂ square. The counterexample check verifies the orthonormal
rows and exact shadow values

| Shadow | Diameter | Smallest square side |
| --- | --- | --- |
| H | 2√(2/3) | 1 + 1/√3 |
| I | 4√2/3 | 4/3 |
| S | 2 | √2 |
| Q | √3 | √2 |

The matrices are given explicitly in the program. H and I already refute the
fixed-coefficient conjecture: Δ = Φ₂₃/₁₀(HO) − (3√2/4)Φ₂₃/₁₀(IO) > 9/1250.
The strict comparison has a rational certificate from
707/500 < √2, 433/250 < √3, √6 < 49/20 and √3 < 7/4.

The square computation includes every orientation. For a centrally symmetric
finite vertex set, the best translation is zero, and the side for axes
(u, Ju) is twice the support function of K = {±v, ±Jv}. On each arc where a
single support vector is active, that positive support function has second
derivative equal to its negative. Its minimum is at an arc endpoint. Those
endpoints satisfy (k − l)·u = 0 for two support vectors. The program enumerates
these exact normals and quadrant endpoints. Thus the calculation is not an
angular mesh. The convex hull does not change the coordinate extrema or the
largest vertex distance.

Additional output checks the algebra for failure at c < c₋ and c > c₊.
It also checks the perturbation budget: assuming the geometric estimates
|diam(A) − diam(B)| ≤ 2d_H(A,B) and |s(A) − s(B)| ≤ 2d_H(A,B),
operator distance at most 1/1000 from H leaves gap greater than 3/5000.
The geometric Lipschitz statements follow by enlarging an enclosing square
by twice the Hausdorff distance and by the triangle inequality; they are
prose arguments here, not Lean results. The preregistered settling target is
the explicit pair H,I. The coefficient intervals and open family belong to
Triage. No assertion is made about sufficiency on [c₋,c₊], or about changing
the source's separately proved Nieuwland number.
