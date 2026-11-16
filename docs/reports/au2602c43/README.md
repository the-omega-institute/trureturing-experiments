# Au coefficient-space Hilbert series: exact finite checks

Source: Kam Cheong Au, *Discovering hypergeometric series with harmonic numbers
via Wilf–Zeilberger seeds*, arXiv:2602.08721v4, Conjecture 4.3 and Example VIII.

Run with Python 3.12 or newer, using only the standard library:

```sh
python3 check.py > check.out
```

`check.py` expands the literal three-parameter seed by its Pochhammer product
recurrence through total parameter degree 8, at every sequence index 0 through
79. All coefficients use `fractions.Fraction`. Gaussian elimination clears
denominators and uses integer gcd reduction, so the ranks are exact over Q.

Both disjoint index windows 0–39 and 40–79 give ranks
1, 2, 4, 6, 9, 12, 16, 20, 25 in degrees 0–8, matching
(floor(N/2)+1)(ceil(N/2)+1). Each finite evaluation rank is a lower bound for
the dimension of the infinite-sequence coefficient space. These checks alone
do not establish the infinite-sequence upper bound or the all-degree conjecture.

Source fidelity checks reproduce the displayed first six Hilbert coefficients
and the exact degree-zero seed identity S(0,0,0;n)=-(1+8n)a_n, where
 a_n=9^(-n)(1/2)_n(1/4)_n(3/4)_n/(n!)^3.
The coefficient substitution (a,b,c)=(x,z,y-z) is computed independently by
binomial expansion and by a transformed product recurrence. The results agree
at all 80 indices, and every coefficient odd in z vanishes through degree 8.

The output is `check.out`; successful execution ends with `ALL_CHECKS_OK`.
The program contains no floating point calculation and no Lean proof.
