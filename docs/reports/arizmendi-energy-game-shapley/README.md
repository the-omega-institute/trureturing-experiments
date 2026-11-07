# Graph energy game: exact Shapley check

Source: Gerardo Arizmendi and Octavio Arizmendi, *The graph energy game*,
arXiv:2012.14114v2, Example 3.9 and Conjecture 3.10;
Discrete Applied Mathematics 330 (2023), 128–140,
DOI [10.1016/j.dam.2023.01.030](https://doi.org/10.1016/j.dam.2023.01.030).
The definitions use ordinary, uncentered adjacency energy.

Run with Python 3 and SymPy (checked with SymPy 1.14.0):

```sh
python3 check.py
```

`check.out` contains the output; success ends with `ALL_CHECKS_OK`.
Every comparison uses integers, rational numbers or exact symbolic algebra.
Radical enclosures use integer square roots at denominator `10**16`, with
both squared endpoints checked. No floating point supports any conclusion.

The program independently reproduces Example 3.9 from the adjacency
matrices of all induced subgraphs of the three-vertex path and the literal
factorial Shapley definition. It checks all eight core inequalities there.
It also compares the complete subset sum with the two orbit formulas on
`K_(2,m)` for `1 <= m <= 4`, checks 30 characteristic polynomials (including
empty parts), and checks both predecessor-probability identities for
`0 <= m <= 70`.

For the family calculation, put `q = sqrt(2)`, `H = (m+1)(m+2)`,
`S1 = sum(sqrt(k), k=0..m)` and `S3 = sum(k*sqrt(k), k=0..m)`.
A scarce-side vertex has Shapley payoff

```
phi_a = 2*((m+q)*S1 + (q-2)*S3)/H.
Delta_m = phi_a - 2*(q-1)*sqrt(m).
```

The analytic input is monotone integral comparison. Comparing each unit
interval with its endpoint values gives

```
(2/3)*m^(3/2) <= S1 <= (2/3)*m^(3/2) + sqrt(m),
(2/5)*m^(5/2) <= S3 <= (2/5)*m^(5/2) + m^(3/2).
```

Since `q-2 < 0`, the first lower bound and second upper bound imply

```
Delta_m/(2*sqrt(m)) >= P(m)/H,
P(m) = ((13-9*q)/15)*m^2 + (1-4*q/3)*m + 2-2*q.
```

The checker verifies the algebra, `q < 99/70` by positive squaring, and

```
P(m) > (19*m^2 - 930*m - 870)/1050,
19*(t+50)^2 - 930*(t+50) - 870 = 19*t^2 + 970*t + 130 > 0 (t>=0).
```

Thus the paper calculation yields a blocking coalition `V \ {a}` for
every integer `m >= 50`. This is an infinite-family certificate using the
stated analytic input, rather than an inference from bounded sampling.
Additional direct rational enclosures certify the deficits at `m=50,64,100`;
the B-vertex payoff is calculated independently and efficiency is checked.
The actual 52-vertex adjacency matrix is checked for simplicity,
connectivity, minimum degree 2, rank 2, `A^3=100*A`, and `tr(A^2)=200`.
The same integral comparisons squeeze the normalized sums, yielding
`Delta_m/sqrt(m) -> (26-18*sqrt(2))/15 > 0`.

These checks are independent exact computational evidence and symbolic
checks of the paper argument. They are not a Lean proof or a claim of
exhaustive literature novelty. The conjecture, its negation, the spectral
bridge on arbitrary finite carriers, Shapley counting and integral
comparisons still require the later formalization stage.
