# LCM splitting and adequacy: exact projection checks

Run `python3 check.py` with Python 3; only the standard library is used. The captured result is in [output.txt](output.txt). All arithmetic uses integers and `fractions.Fraction`.

The fidelity check reproduces the eight signed integer splittings for 6 and 10 and the failed greedy choice for 12 and 18 from Thomas Preu, *LCM Splittings in GCD Domains*, [arXiv:2610.23799v1](https://arxiv.org/abs/2610.23799v1), Section 2. The remaining checks test the proposed converse-refutation mechanism:

- superlevel-set and lattice-operation bounds on rational values;
- disjoint LCM splitting on every ordered pair of three-coordinate vectors over a five-value alphabet;
- sixteen families of bounded functions on positive-integer pairs, with explicit finite support counts for their meets with beta and a diagonal divisor.

The diagonal families use beta(m,n) = 1/(m+s)^p, with s in {0,2} and p in {1,2}. A proposed coprime-factor representative is C(m,n) = 1 for n <= k(m), and 1/(m+n+1) otherwise, for four nonnegative integer profiles k. The selected diagonal is n = k(m)+1, with H = 1/2 there and zero elsewhere. The program checks H <= 1-C, exact finite superlevel certificates for min(C,beta) and min(H,beta), and exactly M entries of H >= 1/2 in M selected rows.

These are finite projections and explicit support certificates for infinite formulas. They neither prove the assertion for arbitrary representatives nor construct a Bézout domain. In particular, on a finite coordinate set every function belongs to the finite-superlevel ideal; the finite quotients are trivial and are not counterexamples. A complete mathematical argument must treat the infinite quotient and the Krull–Kaplansky–Jaffard–Ohm realization theorem. No Lean proof is supplied here.
