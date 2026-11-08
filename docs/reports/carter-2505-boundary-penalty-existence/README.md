# Carter boundary-penalty existence checks

`check.py` checks supporting ingredients for both extensions conjectured after
Proposition 10 in §7.2 of Jack Storror Carter, *Existence of penalised likelihood
estimates and posterior propriety of separable prior distributions for Gaussian
precision matrices*, arXiv:2505.20005v2.

Run `python3 check.py`. The only dependencies are Python's standard library.
All evidence uses integers and `fractions.Fraction`. Logarithms are enclosed by
rational bounds from the positive atanh series, with an explicit geometric tail
bound after 40 terms. `output.txt` is the program's output.

The fidelity check reproduces the determinant, trace and diagonal estimates on
Proposition 4's unit-null-vector path. The proposed extension is checked through
the trace identity on every vector pair of lengths 1–4 with entries −1, 0, 1;
principal-block inequalities; maximal-minor determinant bounds with their
combinatorial exponents; and range/kernel and trace/off-diagonal inequalities.
The random generator has a fixed integer seed.

For the source's logarithmic special case, 20 rank-one sample covariances with
rational observations have an explicitly calculated maximizer at coefficients
strictly above nullity/dimension. The program checks its stationary equations,
negative definite Hessian, and 120 objective comparisons using certified log
intervals. These rational observations lie in Gaussian support; they are not
Gaussian-distributed draws. For S with every entry 1, rescale precision diagonals
by the two nonzero observation coordinates. Writing the rescaled precision as
[[x,z],[z,y]], set w=√(xy) and ρ=z/w. At fixed w,ρ, x+y≥2w gives the maximum at
x=y=w. The remaining objective is
2(1−c)log(w)+log(1−ρ²)−2(1+ρ)w.
Its unique maximum has ρ=−(1−c)/c and
w=c(1−c)/(2c−1), for 1/2<c<1. This justifies calling these analytic maximizers;
finite comparisons alone would not certify a global maximum.

The conjecture uses strict inequalities. At equality c=1/2 it promises no
attainment: [[t+1/4,−t+1/4],[−t+1/4,t+1/4]] has objective increasing to −1.
For c=1/3, the source's path I+t[[1,−1],[−1,1]] has constant trace and
exp(3(F+2))=(1+2t)^3/(1+t)^2≥2t, verifying subthreshold divergence.
Thus a demand for universal attainment at equality would exceed the source.

The final check uses the genuinely mixed family
[[t+1,−t,0],[−t,t+1,0],[0,0,t^(−q)]] with t=2^k and q=1,2,5,
S=[[1,1,0],[1,1,0],[0,0,1]], and a continuous diagonal penalty equal to
(3/4)log(x) below 1 and (2/3)log(x) above 1. This meets the conjectured Gaussian
threshold, and is outside Proposition 10's large-x hypothesis. The checks
certify decreasing objective values and the proposed logarithmic upper bound.
This S is a deterministic stress test of the block reduction, not a claim that
its nongeneric null space represents an almost-sure Gaussian event.

Finite checks do not establish the universal conjecture, compactness, Gaussian
almost-sure general position, or any Lean theorem. The main proposed settlement
is a general analytic argument; these computations are supporting checks only.
No third-party program source was copied.
