# Compact-support Dirichlet spectral gap: exact helper checks

The target is the nonnegative compact-support clause in Kerner and Täufer,
*On the spectral gap of one-dimensional Schrödinger operators on large intervals*,
arXiv:2012.09060v2, Remark 2.10 (journal Remark 10, DOI
10.1007/s00013-024-02060-3). The planned theorem is a general analytic argument
for every fixed nonzero nonnegative bounded measurable compact barrier.

Run with Python 3 and SymPy 1.14.0:

```sh
python3 docs/reports/compact-support-spectral-gap/check.py
```

`output.txt` is its exact output. All evidence uses rational arithmetic or exact
symbolic algebra; no floating-point approximation is used.

The fidelity checks reproduce the differential-expression identity behind
Proposition 2.4's unitary scaling and the source's free Dirichlet gap. The
settlement checks verify the two-branch quantization expansion symbolically,
including an arbitrary energy derivative of the boundary eigenvalue. They then
check every one of 942 piecewise-constant barriers with one to three blocks:
block heights are κ² with κ in {1,2,3}, block widths are log(q)/κ with q in
{2,3}, and gaps between blocks are 0 or 1/2. These are genuine continuum
potentials whose zero-energy transfer matrices have rational entries. There
are 834 nonsymmetric transfer matrices. For every case the program verifies
the transfer determinant, positive boundary response, strictly positive
reciprocal eigenvalue separation, and invariance under equal free padding.
The six one-block cases also recover the symmetric step coefficient.

The predicted coefficient is

\[
C_v=16\pi^2\frac{\sqrt{(A-D)^2+4}}{C},
\qquad
\Gamma_v(L)=C_vL^{-3}+O_v(L^{-4}),
\]

where the zero-energy core transfer matrix is [[A,B],[C,D]]. The exact checks
support identities in this route; they do not prove existence, completeness
or asymptotics of the native infinite-dimensional spectrum, or establish
that the conjecture remains open. Finite positive checks are helpers and do
not supply the planned theorem's admission basis.
