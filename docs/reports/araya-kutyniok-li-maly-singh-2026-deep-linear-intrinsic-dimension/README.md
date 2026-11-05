# Deep-linear random-affine intersection checks

Source: Ernesto Araya, Gitta Kutyniok, Weilin Li, Johannes Maly and Manjot Singh,
*Formalizing the Intrinsic Dimension of Objective Landscapes by a Generalized Notion of
Co-dimension*, arXiv:2610.39674v1, Section 4.1 (equation (20) and the prediction following
Proposition 21), with intrinsic dimension defined in equation (3).

For one hidden layer the parameter set is Θ = {(W,a): Wᵀa=u}. The proposed obstruction
uses q(W,a)=aᵀWe, where e=u/||u||. Its inertia is (M,M,M(D−1)). Near a d-dimensional
subspace of the negative eigenspace, d≤M, q(t)≤−||t||²/4. If ||g||²<||u||/3, then
q(g+t)≤3||g||²/2<||u||/2 throughout that affine plane. Hence the plane misses Θ with
a uniform population-loss gap. The event is open and nonempty under the source's
independent Gaussian-column/Gaussian-offset model. The proposed threshold is at least
max(D,M+1). For D=1, the proposed complete value is M+1, since any (M+1)-plane meets
both signed M-dimensional eigenspaces, giving both unbounded signs on every translate.
These are mathematical outlines, not a formal proof.

Run with Python 3, SymPy 1.14.0 and NumPy 2.4.2:

```sh
python3 check.py
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 monte_carlo.py
```

`check.py` uses exact integers, fractions and symbolic algebra. It reproduces Section
7.2's population-MSE identity on 36 deep networks and symbolically, checks 24 exact
inertia calculations and 216 rational translated-plane maxima, and checks the scalar
square-completion, neighborhood margins and D=1 dimension inequalities. It exits 0
and ends with `ALL_CHECKS_OK`. These finite checks do not establish probabilities or
an infinite theorem. The output is `check-output.txt`.

`monte_carlo.py` samples independent standard Gaussian p-by-d matrices, takes their
Haar column spans by QR, and samples independent standard Gaussian offsets. The
projected offset in Definition 2 generates the same affine set g+U. The teacher is
u=12e₁ throughout; this is a specified member of the source's parameter family.
`monte-carlo-output.json` records counts from 100000 draws per small architecture,
including D=1,2 and several widths, plus 5000 draws at (D,M,d)=(20,30,20). The strict
event is λmax(QᵀBQ)<−1/4 and ||g||²<4. A broader sufficient miss event uses the
quadratic supremum on negative-definite restrictions. These floating-point readings
are diagnostics only, and are separate from the exact arithmetic evidence.

Figure 3 records a numerical-loss crossing at a finite set of sampled planes; it
does not test almost-sure intersection. The large-architecture diagnostic does not
reproduce its unspecified teacher coefficients or initialization: zero observed
misses cannot imply zero probability. The source's tolerance 10⁻⁶ and its first
successful dimension cannot rule out a positive-probability set of failing planes.
The exact proof outline supplies that missing distinction. No Lean proof is included.
