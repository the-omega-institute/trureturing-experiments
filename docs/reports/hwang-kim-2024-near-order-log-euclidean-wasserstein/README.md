# Hwang–Kim: log-Euclidean versus Wasserstein near order

This report checks the final unnumbered open inequality in §6 of Jinmi Hwang and
Sejong Kim, *Near-order relation of power means*, arXiv:2407.07438v1.
The expected answer is negative. These are exact algebraic certificates outside
Lean; the report contains no Lean proof of the assertion or its negation.

Run `python3 check.py` with Python 3 and SymPy 1.14.0. Output is in `check.out`.
Every algebraic constant uses integer/rational arithmetic or exact SymPy radicals;
no floating point is used as evidence. A successful run ends in `ALL_CHECKS_OK`.

The fidelity checks reproduce the source’s congruence formula (2.4) as a
polynomial identity, and its harmonic–LE–arithmetic near-order comparison (6.21)
on a commuting positive diagonal pair. The settlement check uses
A = diag(4,1), B = [[19,−3√15],[−3√15,61]]/64. It verifies positive definiteness,
noncommutation, KAK=B, the binary Wasserstein formula, the spectral calculation
of L=LE, det L=1 and the exact value of √det W. Integer powers certify rational
brackets for u=2^(1/4), z=√265. Rational interval arithmetic then gives
3/2500 < T < 121/100000, where
T=(tr W + d tr L)^2 − (1+d)^2(tr(LW)+2d), d=√det W.

For positive L,W with det L=1, put C=L⁻¹#W. Cayley–Hamilton gives
C=(W+dL⁻¹)/√(tr(LW)+2d), det C=d. The positive T bound implies
tr C>1+d and det(C−I)<0, hence C is not ≥I. This is the exact witness.
The script verifies the matrix polynomial identity behind this reduction.

For each r>1, let x=log r, A=diag(r²,1),
Bθ=Rθ diag(r⁻²,1) Rθᵀ. The script derives the second-order jets from the literal
Riccati equation and the matrix exponential and checks

det(Cθ−I)=D(r)θ²+o(θ²),
D(r)=(-r²+2rx+1)S(r,x)/(32r³(r+1)²),
S=r⁶+2r⁵x−16r⁴x+5r⁴+12r³x−16r²x−5r²+2rx−1.

The analytic completion is as follows. Positive-definite matrix square roots,
inverses and exponentials are real analytic, so the verified jets have the
stated Taylor remainders. Also -r²+2rx+1=-2r(sinh x−x)<0 for x>0.
The identity
S/r³=2sinh(3x)+10sinh x+4x cosh(2x)−32x cosh x+12x
has zero linear coefficient and cubic coefficient 16/3!; the numerator of each
coefficient x^(2k+1)/(2k+1)! for k≥2 is
2·3^(2k+1)+10+4(2k+1)·4^k−32(2k+1)>0.
Indeed 4^k≥16, so the last two terms together are at least 32(2k+1)>0.
Absolute convergence of the hyperbolic series proves S>0 and D(r)<0.
Thus some ε(r) in (0,π/2) makes every 0<|θ|<ε(r) a counterexample, with
strictly positive noncommuting inputs. This all-r analytic argument is prose;
the program verifies its symbolic jets and coefficient identities and also
checks the integer inequalities for k=2,…,128. Those finite checks alone are
not evidence for the all-r statement. The concrete witness already refutes the
source’s universal assertion independently of the family extension.

The separate entrywise-eigenvalue question is not settled. The source’s weak
log-majorization theorem is compatible with this near-order obstruction.
