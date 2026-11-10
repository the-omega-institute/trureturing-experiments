# Exact checks for the all-width extension in Remark 3.12

Source: Eloi Martinet and Leon Bungert, *On the (non-)closedness of the set of ReLU networks*, arXiv:2611.09887v1, submitted 6 November 2026. The target is the extension in Remark 3.12 of Theorem 3.1 to every first width N ≥ 4 and second width M ≥ 6. The program is original code implementing the paper's mathematical formulas and the proposed asymmetric extension; it incorporates no third-party program source.

Run with Python 3.12 or newer:

```sh
python3 check.py
```

The program uses only integers and `fractions.Fraction`. Every polygon vertex, activation decision, determinant, maximum, integral, and comparison is rational. It partitions the full rectangle Q_L = [−L,N+1+L] × [−L,L] by the first-layer and second-layer activation lines of both functions. It verifies the resulting affine formulas against a separate direct evaluator. The maximum of the affine difference is taken at polygon vertices. Squared function error is integrated by triangle barycentric formulas, and squared gradient error by exact polygon areas. Thus these are integrals over the rectangle, rather than quadrature or samples.

Fidelity checks reproduce the source's four-wall reflection identity, the exact active-half coordinate reconstruction, noncollinearity of the four base points, and the uniform estimate 3(λ+2)Lη for the source family, with λ = 16, L = 3 and c_j = j/100. Three source-family convergence instances are integrated exactly.

For the asymmetric candidate, the program checks N = 4,…,12. It enumerates all triples among the 3N designated kink points and all nonvertical lines determined by two points. It verifies that a nonvertical line contains at most two points, and that the coverage inequality forces exactly N designated vertical lines. It checks a nonaffine odd part of the actual target at each designated point and the literal approximating network on a rational grid. Grid observations are finite diagnostics; they do not prove local nonshallowness or global nonrepresentability.

For N = 4, L = 3, η = 1/(96·2^k), k = 0,…,5, exact polygon calculations give:

- uniform error divided by η: exactly 4617/32;
- squared gradient L² error divided by η: strictly between 41000 and 42000;
- squared W¹,² error divided by η: less than 43000;
- maximum squared gradient of the approximants: less than 2400.

The stored output contains the unreduced descriptions as exact rational results and ends with `ALL_CHECKS_OK`. The program checks additional bounds derived from wall-strip and graph-band areas. These finite checks support the general proof outline; they do not establish the assertion for all widths, domains, dimensions, or Sobolev exponents. The planned positive settlement rests on the general incidence and thin-slab arguments.
