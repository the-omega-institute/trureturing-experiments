# Robin cuboid Riesz threshold checks

Source: M. Baur and S. Larson, *Optimizing Riesz means of Robin Laplace operators on cuboids in a semiclassical limit*, arXiv:2604.11076v1, Conjecture 7.1. Source URL: https://arxiv.org/abs/2604.11076v1.

Run `python3 check.py > check.out` with Python 3 and SymPy. The program uses integers, fractions, exact symbolic differentiation, and rational interval arithmetic rounded outwards on a grid of denominator 10^28. Decimal endpoints in the output are rounded outwards too. No binary floating point is used. The last line on success is `ALL_CHECKS_OK`.

Fidelity checks reproduce the exact Weyl-zero value beta_W(1/2,0)=3/4 from Lemma 7.2(ii), using a rational antiderivative after the substitution s=2u/(1+u^2), and verify the strict interval-frequency bounds of Lemma 2.3 in 24 cases. Roots are bracketed by the increasing phase equation l*k+2*atan(k/eta)=pi*j. Atan and cos use alternating-series or Taylor remainder bounds, pi uses Machin's identity, and square roots use integer square roots.

The proposed new estimate is the finite-positive-Robin expansion

R_(alpha,b)(t) = L_(alpha,1) t^(2alpha+1) + B_alpha(b) t^(2alpha)/2 + t^alpha G_alpha(t+2atan(1/b)) + O(t^(alpha-1)),

where G_alpha(x)=Gamma(alpha+1)/pi times the sum over m>=1 of m^(-alpha-1) cos(2mx-pi(alpha+1)/2). The program records normalized remainders and t times their difference from G at three increasing cutoffs. It checks d=1,2 and gamma=1,2: alpha=gamma+(d-1)/2. Integer-alpha Fourier series use exact Bernoulli formulas; half-integer-alpha series have an explicit integral-test tail bound. The asserted finite error bound of 10/t is only a check on these selected inputs, not a uniform theorem.

For d=2 the actual rectangles have sides (1,10^8). The lower and upper rectangle rules and the transverse Dirichlet/Neumann brackets give

abs(R_rectangle/ell - L_(gamma,1) R_(gamma+1/2,b)(t)) <= R_(gamma,b)(t)/ell.

The output includes the resulting enclosure for the finite rectangle's normalized remainder. This checks the dimensional-transfer setting, not a proposed expansion on a fixed square. The longitudinal spectrum is computed from the Robin root equation; the transverse spectrum is enclosed by its Dirichlet/Neumann bounds.

Four further samples exhibit strict Weyl-bound violations at rational b strictly above certified brackets for beta_W, for both dimensions and both gamma values. The rectangle samples include explicit finite-ell error deductions. The root equation, separation of variables, and elementary integral comparison are stated analytic inputs, not kernel-verified code. These checks support the proposed route and certify their finite arithmetic readings; they do not prove the universal conjecture, spectral realization, or the uniform endpoint estimate. There is no Lean proof in this report.
