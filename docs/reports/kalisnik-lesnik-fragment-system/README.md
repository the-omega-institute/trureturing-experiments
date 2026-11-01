# Fragment System exact checks

Source: Sara Kališnik and Davorin Lešnik, *Continuity of Magnitude at Finite
Subsets of ℓ₁ᴺ*, arXiv:2609.26560v1, Question 8.1, Definitions 3.6 and 5.1.
The conventions include dotted fragments from Definition 4.5.

Run with Python 3 and SymPy 1.14:

```sh
python3 check.py > check.out
```

The program uses exact rational and polynomial arithmetic. It reproduces
Example 7.2's two-point formula and Lemma 5.2 at zero on 27 configurations.
It checks the literal dotted-fragment definition against its coordinatewise
implementation, certifies 108 determinants positive by exact Bernstein coefficients after
removing endpoint factors, checks 205,160 instances of Lemma 4.7's distance
identity, and integrates the one-dimensional interval weight measure symbolically
in all three position regimes. The configurations include every nonempty subset
of a three-point one-dimensional grid and a two-by-two planar grid, a skew pair,
a planar cross, and three singleton or pair configurations in dimension three.
Every configuration is tested at four rational radii. The zero-dimensional
empty matrix is checked separately.

For rational distances, a common denominator M converts every exponential
coefficient into a monomial in t = exp(-1/M). The polynomial certificate is valid
for every 0 < t < 1 and uses no floating-point approximation. The output ends in
ALL_CHECKS_OK. These are finite checks supporting a proposed general argument;
they do not settle the universally quantified question or certify Fourier energy
positivity for arbitrary signed measures.
