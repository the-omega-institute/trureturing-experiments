# Positive symmetric unit-step weights and fractional revival

This report checks the existence-classification candidate attached to Duaa Abdullah,
*Quantum Fractional Revival and Entanglement Entropy in Unitary Cayley Graphs*,
arXiv:2605.13645v1, Discussion and Open Problems, item 3.

The candidate answer for order n >= 2 is: a positive real symmetric weight on every
unit step, and zero on every nonunit step, admits two-vertex fractional revival
with nonzero transfer amplitude exactly when n is a power of two or twice an odd
prime power. Perfect state transfer counts under the source's Definition 1.
The full argument and binding definitions belong in the preregistration draft;
these checks are finite mathematical evidence, not a Lean proof or a novelty claim.

Run with Python 3 and SymPy:

```sh
python3 docs/reports/abd-2605-weighted-unitary-cayley-qfr/check.py
```

The stored output is `check.out`. Successful execution ends in `ALL_CHECKS_OK`.
The program uses integer polynomial coordinates in Q(zeta_N), exact rational
bounds and symbolic algebra throughout; it never uses floating-point evidence.

Checks:

- Source fidelity: the unweighted n=8 spectrum forces identical revival phases;
  at n=6 and t=2*pi/3 the transfer probability is exactly 3/4.
- The supplied positive n=8 weights have indexed spectrum
  (4,1,0,-1,-4,-1,0,1), giving the antipodal shift at time pi.
- At n=12 every symmetric positive unit weighting has parameters a,b on the
  steps +/-1 and +/-5. The exact spectrum gives lambda_0=2*lambda_2 and
  lambda_3=0. Equal even phases z satisfy z=z^2; z is invertible, so both
  parity phases are 1. This excludes all positive a,b and all times, rather
  than a finite sample of weights or times.
- For every permitted n<=64, construct an integer spectrum, invert and recompute
  its DFT exactly, check symmetry, reality and vanishing on every nonunit,
  certify strict positivity by a rational bound on the perturbation, check
  every prime-coset sum, and compute the full revival column and normalization.
- For every excluded n<=64, check the necessary phase-order obstruction (or
  the absence of a distinct involutive target for odd n). This tests the
  arithmetic obstruction; it is not an exhaustive search of continuous weights.

For n=2^h, h>=3, perturb the unweighted spectrum on odd indices by the symmetric
anti-half-period sign pattern with entries +/-1. The inverse perturbation has
absolute value at most 1/2 and vanishes on even steps; at time pi the phase is
(-1)^k. Orders 2 and 4 use their ordinary unweighted adjacency.

For n=2*p^a, let m=p^a and L=p^(a-1). Select the L centered residues modulo m,
one from each coset modulo L. Set f(k)=p-1 on those residues and -1 elsewhere.
Each p-coset sum vanishes; f is symmetric. The m-point inverse DFT delta is
real, supported on units, and bounded by 2*(p-1)/p. Add p to every unit weight.
Its integer spectrum p*c_m(k)+f(k) is -1 modulo p. Lifting through
Z/(2m) = Z/2 x Z/m yields two distinct conjugate phases at time 2*pi/p.

No full parametrization of every qualifying weight function is claimed.
The independent literature screen remains limited by inaccessible primary text
for Wang--Feng, DOI 10.1016/j.amc.2023.128006; this report does not resolve that
prior-art question.
