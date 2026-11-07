# Pelikan QSD shape-convergence counterexample checks

Source: Samuel Brevitt and Rainer Klages, *Quasi-stationary and quasi-ergodic
distributions in the Pelikan random map*, arXiv:2607.19084v1, Section 4.2,
equation (52) and its final unnumbered conjecture.

Run `python3 docs/reports/brevitt-klages-qsd-convergence/check.py` from the
experiments repository. Python 3 and SymPy are required; the recorded run used
SymPy 1.14.0. Output is in `check.out`. All evidence is exact rational arithmetic
and symbolic identities. The program uses no floating-point approximations.

The fidelity check reproduces the source's master-equation mass law and its
Section 4.1 rational generating function for a QSD with initial coordinate
alpha and eigenvalue 1-p*epsilon*alpha. This is independent of the proposed
counterexample coefficients. It uses equation (52) and the master equation;
the displayed evolution equation labelled eq:proevo in the source TeX omits
the q term that is present in those formulas.

Set p=3/4, epsilon=1/2, r=4/5, lambda=73/80. For every k>=0 define

- u_k = (7/6)(1/2)^k + (28/345)(4/5)^k - (70/69)(5/12)^k;
- v_k = -(153/286)(1/2)^k + (18972/1495)(-4/5)^k
  - (1785/253)(-5/12)^k.

These are full infinite sequences, represented by finite sums of geometric
terms. Their rational generating functions are

- U(z) = (7/30)(1-3z/4)/[(1-z/2)(1-4z/5)(1-5z/12)];
- V(z) = (51/10)(1-3z/4)/[(1-z/2)(1+4z/5)(1+5z/12)].

All geometric ratios have absolute value less than one, so summing them as
c/(1-s) is valid. The sums of u and v are both one. Here v is a signed
sequence, not a probability distribution. The linear transfer T, which is
the numerator of equation (52), satisfies Tu=lambda*u and Tv=-lambda*v.
The program verifies the boundary coordinate and each geometric-mode
coefficient, and independently verifies the whole rational-function identity.
It also evaluates coordinates 0 through 100 as a separate arithmetic check;
this finite check does not establish the all-index result.

The inequalities (1/2)^k >= (5/12)^k and |s|<=r for all appearing ratios give
u_k >= (28/345)r^k, u_k<=3r^k, and |v_k|<=21r^k. For every real delta in
0<delta<=1/1000 define a_delta=(u+delta*v)/(1+delta) and
b_delta=(u-delta*v)/(1-delta). Their numerator lower bound is
(4151/69000)r^k>0, and both normalized sequences are bounded above by 4r^k.
Thus a_delta, b_delta and u are strictly positive probability distributions.
Their k-th roots converge to r by squeezing between positive constant
multiples of r^k. This is a mathematical implication of the certified bounds,
not a numerical root-limit test. All three have shape sigma=log(5/4).

The normalized operator F has survival denominator 1-(3/8)x_0. The exact
symbolic identities give F(u)=u, F(a_delta)=b_delta and F(b_delta)=a_delta.
Their survival masses are lambda*(1-delta)/(1+delta) and
lambda*(1+delta)/(1-delta), respectively. Their coordinate-zero gap is
146*delta/[15*(1-delta^2)]>0. Consequently the even and odd subsequences of
the iterates are distinct constants. The orbit fails to converge pointwise,
and therefore fails in l1, total variation and weak convergence on the
discrete state space. This implication is part of the written argument;
no Lean proof or kernel-checked settlement is claimed by this experiment.

At delta=1/1000, a_0=7153/30030, b_0=6847/29970, the gap is
29200/2999997, and the survival masses are 72927/80080 and 73073/79920.
The raw chain is irreducible because both adjacent transitions are positive,
and aperiodic because its self-loop at zero has probability 3/16. Its period
is distinct from the two-cycle of the normalized distribution dynamics.
