# Attracting rank cycle for Aldous–Cruz–Feng Conjecture 5

Source: David J. Aldous, Madelyn Cruz and Shi Feng, *Markov chains and mappings
of distributions on compact spaces II: Numerics and Conjectures*,
[arXiv:2403.18153v1](https://arxiv.org/abs/2403.18153v1), §5.1, Conjecture 5.

Run `python3 check.py`. The saved standard output is `check.out`.
Python 3 standard library only; all evidential calculations use integers and
`fractions.Fraction`. Successful completion prints `ALL_CHECKS_OK`.

The fidelity check reproduces the exact four-state rank matrix example R2,
its transition matrix, and its invariant distribution from the source.
An independent enumeration of sampled tuples checks the binomial-tail kernel
on small instances, including repeated samples and distributions with zeros.

The proposed refutation uses the line metric with points 0, 4, 7, 12 and,
more generally, every gap triple satisfying 0 < b < a < c < a+b. The third
nearest state forms a single four-cycle. The program checks the entire cone
by positive linear-coefficient certificates, the exact stationary inverse
identity, and all rational constants in the uniform contraction argument
for (j,k) = (5m,8m), m ≥ 128. It also performs nine exact stationary solves
at large k as independent instance checks. The proposed basin has positive
three-dimensional volume and all its limiting coordinates are at least 15/64.

The uniform attraction conclusion additionally uses the analytical binomial
Hoeffding inequality, its application to the binomial derivative, and the
contraction mapping theorem. Exact computations check the algebra and
constants of that argument; they are not a substitute for these analytical
steps or a kernel-checked proof. No Lean proof of the conjecture or its
negation is included in this report.
