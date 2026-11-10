# Adjacent-dimension smoothing sharpness: exact family checks

Source: Soichiro Suzuki, *An order-interpolation inequality for Bessel functions*,
arXiv:2607.00109v1, Section 3, inequality (3.4) and the following open sharpness question.

Run `python3 check.py`. The program uses integers and `fractions.Fraction` only.
The stored output is `check.out`; successful execution ends with `ALL_CHECKS_OK`.

The fidelity check certifies 27 strict instances of the source's main Bessel
interpolation inequality, using alternating-series enclosures and a rational
Machin-formula enclosure of π. It also checks the angular-order inclusion behind
the source's comparison between dimensions d and d+2.

The proposed sharpness family uses the nonnegative Fejér polynomial F_M(2r) in
even dimensions and F_M(2r−π/2) in odd dimensions. The program checks its exact
sum-of-squares expansion, mean, harmonic budget, and target mean 2−1/(M+1) for
M=1,…,32. It verifies the needed angular index and phase for d=2,…,100, and 153
Taylor coefficients of the Legendre sine representation used for uniform
half-integer-order estimates.

These finite identities support the proposed general analytic argument. They do
not certify the uniform transform estimates, logarithmic asymptotics,
essential-supremum passage, Fourier-core density, or the full sharpness theorem.
No numerical quadrature or floating-point result is offered as evidence.

Classical mathematical inputs: Bessel defining series and asymptotics (DLMF
10.2.2, 10.17.3), the spherical addition theorem (DLMF 10.60.2), and Machin's
identity π = 16 arctan(1/5) − 4 arctan(1/239). The program is independently written;
no third-party program or Lean source is copied.
