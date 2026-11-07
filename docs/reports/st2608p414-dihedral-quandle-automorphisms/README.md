# Integral dihedral quandle algebra automorphisms

`check.py` is an independent exact-integer check of the proposed group description for Sangare–Ta, arXiv:2608.01550v2, Problem 4.14. It uses only the Python 3 standard library; `check.out` contains its complete output.

Run `python3 docs/reports/st2608p414-dihedral-quandle-automorphisms/check.py` from the repository root.

The fidelity check recovers the eight permutation automorphisms and the explicit generators and relations in Bardakov–Passi–Singh's computation of the integral algebra of the four-element dihedral quandle. For m = 1 through 8, the checks verify signed affine matrices and their inverses, all affine parameters in the reduction kernel, its four graph-kernel lifts, integral shear lifts and inverses, all basis products for each tested lift and mixed product, the semidirect multiplication formula, the explicit polynomial-order block equations (including equivalence checks on forty arbitrary blocks per m), and the trace-square identity. The arithmetic uses integers, including exact-division assertions in Bareiss elimination and integral gluing. The tested shear parameters are −2, −1, 1, 3. Mixed-product sampling has fixed seeds.

These checks validate finite instances of the proposed formulas. They do not enumerate the infinite kernel, prove completeness of the automorphism description, establish the universal theorem, or constitute Lean verification. The conjectured universal argument is a preregistration target.
