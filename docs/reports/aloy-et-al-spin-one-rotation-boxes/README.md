# Spin-one rotation boxes: exact checks

Source: A. Aloy, T. D. Galley, C. L. Jones, S. L. Ludescher and M. P. Müller,
*Spin-Bounded Correlations: Rotation Boxes Within and Beyond Quantum Theory*,
arXiv:2312.09278v2; Communications in Mathematical Physics 405, 292,
DOI: 10.1007/s00220-024-05123-2. The target is the unnumbered question at the
end of §4.2.3: whether Q₁ᴬ = R₁ᴬ for finite outcome sets of size at least three.

Run with Python 3 and SymPy (checked with Python 3.14 and SymPy 1.14.0):

```sh
python3 check.py
```

The program uses integers, rational numbers, and exact symbolic expressions.
Its output is stored in `output.txt`.

- Fidelity: reproduce the explicit quantum realization of sin⁴(θ/2) in
  Appendix C.3, proof of Theorem 6, item 2, including positivity and POVM completion.
- Symbolic identity: for an arbitrary Hermitian 3 × 3 matrix Q, reversal-transpose
  averaging B = (Q + J Qᵀ J)/2 preserves φ*Qφ at the same angle, where
  φ = (e⁻ⁱᶿ, 1, eⁱᶿ). It forces B₀₁ = B₁₂ and B₀₀ = B₂₂. Fourier normalization
  then forces the sum of the outcome matrices to be diagonal.
- Singular support: check the whitening, POVM-sum, and Born identities on all
  seven nonempty supports, with arbitrary symbolic positive diagonal entries.
- Finite catalog: enumerate every multiset of three or four rank-one matrices
  vv*, with v in {0, ±1, ±i}³ modulo global unit phase. Zero factors and repeated
  factors are allowed. Of 5,984 three-factor and 52,360 four-factor multisets,
  respectively 91 and 575 have constant nonzero sum polynomial. Divide each
  family by its total trace. Check diagonal normalization, support containment,
  principal-minor positivity certificates for the completed effects, and the
  entrywise POVM sum. The symbolic Born identity applies to each such family.

The finite catalog does not exhaust arbitrary positive semidefinite matrices
or nonnegative real trigonometric polynomials. The checks do not prove scalar
Fejér–Riesz factorization or the universal claim in Lean. This report supplies
exact checks for a preregistration draft before any formal proof probe.
