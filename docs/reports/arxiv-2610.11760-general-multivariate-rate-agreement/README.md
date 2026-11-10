# Full-dimensional CDF loss and Wasserstein risk

Source: Catia Scricciolo, *Adaptive minimax multivariate L1-deconvolution:
harmonic-mean rates under noise filtering, and a 1-Wasserstein equivalence*,
[arXiv:2610.11760v1](https://arxiv.org/abs/2610.11760v1), Definition 2.1 and
Section 8. The source asks whether its minimax-rate agreement extends to
general multivariate measures. This report checks the unrestricted family
interpretation with the printed, unweighted Lebesgue CDF loss on all of R².

Run `python3 check.py` with Python >= 3.12 and SymPy 1.14.0. The output is
`check.out`, ending in `ALL_CHECKS_OK`. All certificates use integers,
rationals, exact radicals, symbolic exponentials and symbolic integration.
There is no floating point evidence. The analytic deductions below are
mathematical arguments; this report does not contain a Lean settlement proof.

## Fidelity and admissible family

The program independently integrates the Laplace density to reproduce the
source's characteristic function 1/(1+t²), its reciprocal derivatives 2t and
2, and its variance 2. It also evaluates the source's rate formula at
d=2, alpha=(1,1), beta=(2,2): E=2/5 and psi_n=n^(-1/5).

Let U,V be independent standard normals, theta in [0,1], and

    X_theta = (theta+U, U/2+sqrt(3)V/2).
    Sigma = [[1,1/2],[1/2,1]].

Each observation is X_theta plus an independent vector of independent
standard Laplace variables. Replicates are independent. These signals are
dependent, smooth, have full support, and their covariance eigenvalues are
1/2 and 3/2. Translation preserves all their regularity bounds.

The Gaussian score in coordinate j has second moment (Sigma^-1)_jj=4/3.
Cauchy-Schwarz gives the L1 derivative bound sqrt(4/3)<2, and integration
along a coordinate segment gives ||f_theta(.+h e_j)-f_theta||_1<=2|h|.
For alpha_j=1 the source's m_j=ceil(alpha_j)-1 is zero. Normalization and
this finite-difference bound give N_1²((1,1),2) membership.

For p=1, max(1,X1²)max(1,X2²) <= (1+X1²)(1+X2²), whose exact expectation
is 9/2+2 theta²<=13/2. Thus the source's moment condition holds uniformly.
The Laplace ordinary-smooth bounds have beta_j=2, c0=1, c1=2 and
reciprocal derivative constants D1=D2=2. Both coordinates are active.

## Exact infinite CDF risk

Write F_theta for the signal CDF. For every theta<eta and every real x,

    lim_(y -> infinity) (F_theta(x,y)-F_eta(x,y))
      = Phi(x-theta)-Phi(x-eta)>0.

Tonelli gives ||F_theta-F_eta||_1=infinity. The program additionally checks
an explicit strip minorant for theta=0, eta=1. For x in [0,1/2], y>=1/4,
the rectangle Z1 in [-1/4,0], Z2 in [0,1/4] is contained in
{x-1<Z1<=x, Z2<=y}. The centered Gaussian density is

    (pi*sqrt(3))^-1 exp(-(2/3)(z1²-z1*z2+z2²)).

Its value on that rectangle is at least exp(-1/8)/(pi*sqrt(3)). Therefore
F0-F1>=exp(-1/8)/(16*pi*sqrt(3)) on the strip. Integrating over strips
of arbitrary height T proves divergence, using a checked positive slope.

Every finite-sample observation law has a strictly positive Lebesgue
density, so the laws at theta=0 and theta=1 are equivalent. For any
measurable CDF estimator, either its risk at 0 is infinite, or its loss
to F0 is finite almost surely at 0 and hence at 1. Infinite pairwise
separation and the extended triangle inequality then force infinite loss
to F1 almost surely at 1. Consequently R_CDF(n)=infinity for every n>=1.
The same argument permits arbitrary additional independent randomization:
equivalent observation laws remain equivalent after the same Markov kernel.
The equality includes the upper bound, since infinity is the largest
extended nonnegative value; there is no finite CDF minimax rate here.

## Matching Wasserstein bounds

Translation coupling gives W1(mu_theta,mu_eta)<=|theta-eta|; the first
coordinate mean gives the reverse bound for every coupling. Hence equality.
Clip the first-coordinate observation mean to [0,1] and return the signal
law indexed by that estimate. Clipping cannot increase parameter error;
Var(Y1)=1+2=3 gives the uniform upper bound sqrt(3/n).

For the lower bound take theta0=0 and theta1=delta=1/(2*sqrt(n)). Both lie
in [0,1] for n>=1. Gaussian sample KL is

    n delta² (Sigma^-1)_11 / 2 = 1/6.

Convolution with the same Laplace noise decreases KL. Pinsker therefore
bounds observation total variation by sqrt(1/12)<1/2. The nearest-law test
obtained from any W1 estimator and the metric triangle inequality give
max risk >= delta(1-TV)/4 >= delta/8=1/(16*sqrt(n)). This also covers
randomized estimators. Thus for every n>=1,

    1/(16*sqrt(n)) <= R_W1(n) <= sqrt(3/n).

These bounds determine the rate n^(-1/2) on this actual dependent family;
they are not a comparison between two unmatched bounds.

## Scope

The strip obstruction also applies to products with different marginals.
Thus the literal norm challenges the printed product-class premise too.
This report neither endorses that premise nor computes minimax rates on
the entire unrestricted Nikol'skii class. It does not address losses on
a bounded observation domain, weighted CDF losses, fixed-marginal classes,
or the sum of marginal one-dimensional CDF losses. Those are different
questions, and their agreement is not refuted by this argument.
