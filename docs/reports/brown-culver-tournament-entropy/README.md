# Fixed-order tournament entropy checks

This report checks the integer-moment route to refuting Conjecture 1 of David E. Brown, Eric Culver, Bryce Frederickson, Sidney Tate and Brent J. Thomas, *Entropy of Tournament Digraphs*, [arXiv:1812.09458v2](https://arxiv.org/abs/1812.09458v2), Section 2. The source's exponent is a fixed positive integer greater than one. Its normalization is the directed outdegree Laplacian divided by its trace.

Run with Python 3, without dependencies:

```sh
python3 docs/reports/brown-culver-tournament-entropy/check.py
```

The stored output is `check.out`. Every operation used as evidence is an integer operation or a `fractions.Fraction` operation. The program exits zero and ends with `ALL_CHECKS_OK` exactly when all assertions succeed.

The program enumerates every labelled tournament on 2 through 6 vertices. It reproduces the source's Section 2 table: score-sequence counts 1, 2, 4, 9, 22 and distinct quadratic-entropy counts 1, 2, 3, 6, 9. It also checks the source's parity-dependent formula for the latter. Quadratic entropy has the same distinct-value count as the positive integer trace of the squared Laplacian, because normalization is constant at each order and the logarithm is injective. A separate cyclic-product calculation checks traces through exponent four at orders at most four.

For exponents 2 through 8, the enumeration checks the integer trace bound `abs(trace(L^k)) <= n * (2*(n-1))^k` and the corresponding row-sum norm bound. At order six it finds 36 fourth-power traces, exceeding the 22 score sequences. Thus the source's motivating finite-order phenomenon is compatible with the proposed asymptotic refutation.

The family check constructs actual tournament adjacency matrices from zero through eight three-vertex blocks and zero, one or two terminal singleton blocks. Each three-vertex block independently uses a triangle or a transitive orientation. Later blocks beat earlier blocks. It computes the resulting outdegrees, checks the predicted sorted score lists, recovers every input bit, and checks exactly `2^m` distinct lists for every tested order. All 1,533 constructions pass.

The positivity-bound check evaluates the proposed rational lower bound at `n = 8^k + 1` for exponents 2 through 12. The tail-majorant check verifies that its value at order 900 is less than `1/2^64` for the same exponents. These are finite checks of formulas, not an empirical proof of an infinite limit. The proposed general argument combines a polynomial bound on the integer moment's value range with `2^floor(n/3)` realized score sequences.

The domain-edge check records the triangle's unnormalized traces 3, 0, -9, -27, -54 for exponents 2 through 6. A real logarithm is undefined at the nonpositive values; the enumeration reports positive-trace counts separately. The proposed proof must establish positivity for every tournament eventually, at each fixed exponent. No Lean proof of the conjecture or its negation is included in this report.

The check program is independently written for this report. Its only external inputs are the cited definitions, table, and counting formula; it embeds no downloaded code or data.
