# Kubo tropical coordinates: inverse stability checks

Source: Susumu Kubo, *Stable complete coordinates for multisets of points via basic r-symmetric tropical polynomials*, arXiv:2606.30184v1, equation (1), Section 3, and Section 4's parity-box example.

Run with Python 3:

```sh
python3 docs/reports/kubo-tropical-coordinates-inverse-stability/check.py
```

The program uses integers and `fractions.Fraction` throughout. Frobenius and Euclidean norm inequalities are compared after squaring nonnegative quantities. It independently evaluates the source's disjoint-subset definition on the parity-box example and reproduces its unequal degree-four basic values, equal degree-three values, and equal pairwise projections.

It checks the proposed bound d(M,N)^2 ≤ 64 n^3 r max_c |b_c(M)−b_c(N)|^2 on every unordered pair of multisets in the printed finite domains. It also checks the source's forward bound, its r=1 inverse bound, the transform identity at rational grid points, finite diagonal slope identities away from coordinates, orthant and cell counts, and a single matching assembled across cells. Rational perturbations of separated clouds exercise negative values, repetitions, zero error, and nontrivial gap cuts. The exact cut-boundary case is included.

`output.txt` records 23,461 pairs, 6,969 transform values and 24,062 gap grids, ending in `ALL_CHECKS_OK`. These finite checks support a proposed general mathematical argument; they do not prove its unrestricted real-domain conclusion or the optimal growth of distortion.
