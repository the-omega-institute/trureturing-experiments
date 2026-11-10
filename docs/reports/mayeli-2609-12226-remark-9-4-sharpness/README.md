# Spatial sharpness check

Run `python3 check.py`. The saved output is `check.out`.

The program uses integer arithmetic and rational fractions. It reproduces the
finite trace/correlation formula in Mayeli, arXiv:2609.12226v1, Proposition 4.1,
and the half-band interval specialization of Proposition 9.3. The normalization
is normalized Haar measure on the unit torus and the kernel
`K(k) = sin(pi*k/2)/(pi*k)`, with `K(0) = 1/2`.

It checks the proposed alternating-interval spatial construction at 30 exact
rational scales and parameters, including endpoint exclusions, separated runs,
translation counts, the proposed trace lower constant and transverse tensor
scaling. Rational bounds for pi come from the alternating arctangent series and
Machin's identity. These are finite checks, not a proof of the universal
asymptotic claim or a claim of literature novelty.

Source: Azita Mayeli, *Trace-defect bounds for discrete Fourier concentration
operators*, arXiv:2609.12226v1, Remark 9.4 (spatially rough off-critical branch).
The program is original code under this repository's licence; no external
program or Lean source was copied.
