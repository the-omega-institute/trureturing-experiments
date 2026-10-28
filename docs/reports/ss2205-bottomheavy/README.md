# Sagan–Swanson coefficient bottom heaviness

`check.py` uses Python 3 and its standard library, with integer coefficients only.
Run `python3 check.py`; `output.txt` is the successful output. The source is
Bruce E. Sagan and Joshua P. Swanson, *q-Stirling numbers in type B*,
[arXiv:2205.14078v1](https://arxiv.org/abs/2205.14078v1), Conjecture 7.6 and
Section 7.5. The targeted polynomials are the unshifted, unsigned families.

The checker reproduces the complete printed coefficient lists of `S_B[6,4]`
and `c_B[7,5]`, verifies their failure of unimodality, and reproduces
`S[4,3] = 3 + 2q + q^2`. It checks the conjectured inequality for 11,628
polynomial cases (all four families, `0 <= n <= 50`, `-3 <= k <= 53`),
including zero and constant polynomials, using 1,187,487 reflected coefficient
comparisons. Coefficient arrays omit trailing zeros and retain exponent zero.

An independent monomial-enumeration implementation checks the four elementary
and complete homogeneous symmetric-polynomial descriptions for `0 <= n <= 8`
and `0 <= k <= n+1` (216 identities). General convolution is cross-checked
against the recurrence implementation's sliding-sum multiplication. Further
finite checks cover 1,555 products of up to four q-integers of lengths 1–6,
161,736 comparisons at larger ambient degrees, and 7,056 nonnegative weighted
sums of two short products.

These are finite checks supporting a preregistration. They do not establish the
universal conjecture and do not constitute a Lean proof. A failed assertion
exits nonzero; success ends with `ALL_CHECKS_OK`. The recorded run used
CPython 3.14.8 on macOS arm64; other platforms were not tested.
