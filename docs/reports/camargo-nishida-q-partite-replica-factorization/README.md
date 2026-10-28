# Higher-party bosonic Gaussian replica factorization

`check.py` supplies exact finite evidence for the proposed affirmative answer to
Camargo–Nishida, *Genuine Multi-Entropy of Fully Symmetric Gaussian States*,
arXiv:2609.30754v2, Section 6, fourth future-directions bullet.

Run with Python 3.10 or newer and SymPy (`python3 -m pip install sympy`):

```sh
python3 check.py
```

The stored output uses SymPy 1.14.0 and ends with `ALL_CHECKS_OK` (exit 0).
There is no floating point arithmetic.

The program constructs the full replica quadratic matrix directly from the
hypercubic twists, and compares its normalized determinant with a product of
purities computed independently from reduced position–momentum covariance
matrices. It also checks the partial-parity determinant for each factor,
complementary-subsystem equality, positivity of all determinants, binary
character orthogonality and translation eigenvalues for replica ranks 1–6,
and displacement phase coefficients for party counts 2–6.

The declared finite family is exhausted: q = 2,3,4,5; N = q or q+1; every
ordered composition of N into q positive party sizes; and the two explicitly
defined real and complex Gaussian precision matrices for each N. All 36 cases
pass. The eight q=3 cases reproduce the Section 4 source result
Z₂⁽³⁾ = E_A E_B E_C. The q=4 and q=5 cases test seven and fifteen factors.

The determinant comparisons check squared partition functions and squared
purities, with positive determinants. The positive Gaussian-integral branch,
normalization, partial-trace bridge and all-real-parameter/all-q statement
remain analytic proof obligations; the preregistration outline supplies their
proposed argument. These finite tests are neither a universal proof nor a Lean
proof. First moments are tested through their permutation/phase coefficients.

Source: <https://arxiv.org/html/2609.30754v2>.
Standard hypercubic twists: Akella–Iizuka–Miyata, arXiv:2607.06050v2, Appendix A,
<https://arxiv.org/html/2607.06050v2>.
