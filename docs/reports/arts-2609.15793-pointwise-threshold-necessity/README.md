# Pointwise threshold necessity in the degree-corrected Poisson SBM

The source is Lucie Arts, *A Threshold Constant for Exact Community Recovery in
the Degree-Corrected Poisson Stochastic Block Model*, arXiv:2609.15793v1,
Section 3.2, Remark 1. The proposed settlement concerns its pointwise
no-algorithm assertion, rather than a minimax assertion over nuisance weights.
The source has different argument orders for H in its main definition of C and
in its appendix. Both orders satisfy the upper bound below.

Run `python3 check.py` with Python >=3.12 and SymPy 1.14.0. The exact certificates
use no floating point. Their output is `check.out`, ending in `ALL_CHECKS_OK`.
The check reproduces the source's conditional block-count and degree-sum
expectations, including doubled self-loops, for all 256 pairs of binary
labelings at n=4. It checks the weight construction on all 2047 binary
labelings at n=0,...,10 and checks the algebraic family certificates
symbolically for L=64+u, u>=0. These are certificates for the algebra and the
finite fidelity checks; a Lean proof of the probability argument is not supplied.

For each fixed integer L>=64, use pi=(1/2,1/2),
S=((L,2L),(2L,4L+1)), and rho_n=log(n)/n eventually. A positive finite prefix
may be chosen as min(log(n)/n,1/[2(4L+1)]) for n>=2 and 1/[2(4L+1)] for
n=0,1, so rho_n S has entries in [0,1] at every n. Within each true block,
assign weights 99/100 and 101/100 to successive vertex pairs and weight 1 to
an unmatched vertex. More generally the degree argument works for every
normalized weight assignment in [99/100,101/100]. The two blocks have the same
permitted weight values. Each fixed matrix is strictly positive with determinant L.

## Analytic argument checked by the algebraic certificates

Write alpha=3L/2, beta=(6L+1)/2, p=(1/3,2/3), and
q=(2L/(6L+1),(4L+1)/(6L+1)). These are the masses and normalized pi-weighted
columns. Let KL_s(x||y)=x log(x/y)+y-x and K(p||q)=sum p_a log(p_a/q_a).
For positive columns u,v of masses alpha,beta,

    sum pi_a KL_s(u_a||v_a) = KL_s(alpha||beta)+alpha K(p||q).

The penalized Chernoff function is concave, equals zero at t=0, and is bounded
above by t times its derivative there. With the appendix argument order,
choose the ordered pair with the lower-mass column in the first argument of H.
Its derivative is alpha K(p||q). With the main-text order, choose the reverse
ordered block pair; the derivative is

    alpha K(p||q)+KL_s(alpha||beta)-KL_s(beta||alpha).

The last difference is nonpositive for beta>=alpha: put r=beta/alpha and use
log r>=2(r-1)/(r+1). The latter inequality follows by differentiating the
difference: its derivative is (r-1)^2/[r(r+1)^2]>=0, with equality at r=1.
The elementary bound log x<=x-1 gives

    alpha K(p||q) <= alpha sum (p_a-q_a)^2/q_a = 1/[12(4L+1)] < 1.

Both definitions of C are nonnegative because t=0 is allowed. The appendix
definition is strictly positive: its derivatives at zero are positive for
both ordered pairs, since p differs from q. No assertion depends on assigning
a zero value to an undefined logarithm or on a singular matrix.

Let E={|n_1/n-1/2|<=1/100}. Since n_1 is Binomial(n,1/2),
P(E complement)<=2500/n. Conditional on the entire labeling, the off-diagonal
degree D_i is Poisson with mean

    mu_i=rho_n theta_i (S_{z_i,1} n_1+S_{z_i,2} n_2-theta_i S_{z_i,z_i}).

On E, for all n>=400 in the logarithmic regime, the mean is at most
(25L/16)log n for class 1 and at least (47L/16)log n for class 2. For the
lower bound, discard the positive additive term 99/100*49/100 and bound the
diagonal correction using theta_i<=101/100; the exact remaining slack at
L=64, n=400 is positive. The checker certifies the polynomial slack for
every L>=64. Classify using only whether D_i<(35L/16)log n.

For a Poisson D with mean mu, the exponential series gives
E[r^D]=exp(mu(r-1)), for either r=5/4 or r=3/4. Markov's inequality gives
the upper-tail exponent at most

    (25L/64)-(35L/16)log(5/4) <= -3L/64,

and the lower-tail exponent at most

    -(47L/64)+(35L/16)log(4/3) <= -37L/384.

Here log(5/4)>=1/5 follows by integrating 1/(1+x), and
log(4/3)<=7/24 follows by integrating
1/(1-x)<=1+x/(1-y) on 0<=x<=y=1/4. A union bound, without independence
between different vertex degrees, yields

    P(any error)<=2500/n+n^(1-3L/64)+n^(1-37L/384) -> 0.

At L=64 the matrix is ((64,128),(128,257)), the threshold is 140 log n,
C<=1/3084, and the two powers are n^-2 and n^(-31/6). The degree rule is a
graph-only estimator with a fixed numerical threshold coefficient. It receives
neither true labels nor degree weights. Exact raw labels imply exact recovery
up to label permutation, as required by the source.

The mechanism is a mismatch between normalized profiles, which are close, and
absolute mean degrees, which remain separated. The bound on C tends to zero
as L increases, while the recovery exponents grow. This does not identify the
threshold of profile likelihood or settle a uniform minimax nuisance problem.

## Simulation diagnostic

Run `python3 sanity.py` with NumPy 2.4.2. The fixed-seed output is `sanity.out`.
It samples actual shared upper-triangular Poisson edges for five graphs at
L=64, n=4096; the two endpoints of an edge receive the same sampled count.
It uses floating point solely as a sanity diagnostic. Its results supply no
proof, exact bound, asymptotic evidence, or pass condition for `check.py`.
