# Petz-iterate norm convergence: exact block checks

Source: Saptak Bhattacharya, *Recoverable states on von-Neumann algebras*,
arXiv:2605.08829v1, equation (2), Theorem 1, Corollaries 3 and 4, and Section 2,
Remarks item 2 after Theorem 2. The question asks whether the Petz iterates converge
in operator norm on L1. This directory contains an exact-arithmetic experiment,
not a Lean proof or an infinite-dimensional certification.

Run `python3 check.py`. Python 3's standard library is sufficient. All evidence
uses integers and `fractions.Fraction`; the stored output ends in `ALL_CHECKS_OK`.

## Model and exact identities

Let j range over the nonnegative integers, let w_j=2^(-j-1), and give each of the
two coordinates in block j trace weight w_j/2. On the bounded sequences define
phi_j=(1-epsilon_j)I+epsilon_j S, where S swaps the two coordinates and
epsilon_j=1/(j+3). The reference state is B=1. Each block is stochastic and
symmetric; thus its trace adjoint equals itself, and equation (2) gives R=phi.
Set T=phi^2 and let E average each two-coordinate block. Then

    lambda_j=(1-2/(j+3))^2,
    T_j^n-E_j = lambda_j^n (I-E_j),
    ||T_j^n-E_j||_(1->1)=lambda_j^n.

The script checks 2,673 residual identities, the normalized antisymmetric
witnesses, the geometric trace weights, and the exact norms on finite initial
block sets. It reproduces Theorem 1 on a finite two-point algebra with mixing
parameter 1/4: both the L1 and L2 residual norms are (1/4)^n. It also checks
Corollary 3's L1 contractivity on 169 rational inputs.

For the first N blocks, the exact operator norm is

    ||T_N^n-E_N||_(1->1)=(1-2/(N+2))^(2n).

For fixed N this tends to zero. To test the infinite-family mechanism, take
N(n)=8n+1. The last included block has j=8n, and Bernoulli's inequality gives

    (1-2/(8n+3))^(2n) >= 1-4n/(8n+3) > 1/2.

The script checks these inequalities for n=1,...,128, as well as explicit
near-1 witnesses. This is a moving truncation; it does not assert a positive
lower bound for all times at one fixed finite N.

## Infinite argument and limits of the checks

The diagonal representation on l2 realizes the bounded sequences as a
commutative von Neumann algebra. The positive summable weights give a faithful
normal tracial state. Each block map is a convex combination of the identity
and a swap automorphism; the coefficient is central and varies by block.
Hence phi is completely positive, unital, trace preserving, and strict.
The fixed points of T are exactly the block constants, and E is their
trace-orthogonal projection.

For each fixed n, lambda_j^n tends to 1 as j tends to infinity. The normalized
antisymmetric block vectors give the lower bound 1, and the block formula gives
the upper bound 1, on the infinite L1 operator norm. Every fixed L1 input still
converges: truncate its summable weighted absolute-value tail, converge on the
remaining finite blocks, and use ||I-E||_(1->1)=1 to control the tail.
These limiting and operator-algebraic steps are mathematical arguments to be
formalized after preregistration; they are not certified by finite checks.

The same family also prevents L2 operator-norm convergence. In Theorem 1's
proof, the inference from positivity, contraction, and a nonzero fixed space to
isolation of 1 in the spectrum is invalid without an additional hypothesis.
The fidelity calculation above concerns the valid finite-dimensional instance
of that theorem and does not endorse its general assertion.
