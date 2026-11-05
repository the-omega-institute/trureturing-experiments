# Gauss–Kuzmin permutation-collision density

Exact finite checks for Shreyas Singh, Zhuo Zhang and AJ Hildebrand,
*An Elementary Characterization of the Gauss–Kuzmin Measure in the Theory of
Continued Fractions*, arXiv:2511.01992v1, Integers 25 (2025), A89, Conjecture 6.

Run `python3 check.py` with Python 3. The program uses the standard library,
integer polynomial coefficients and `fractions.Fraction`. `output.txt` contains
its output. Successful execution ends with `ALL_CHECKS_OK` and exit status zero.

The source matrix product and an independent scalar continuant recurrence agree
on all 5,460 words of lengths 1 through 6 over digits 1 through 4. The determinant,
reversal and exact endpoint-ratio formulas hold on that range. The program
reproduces Tables 4–6 at (N,n) = (10,4), (20,4), (10,5), (10,6), including
f = 10, 30, 8, 23. It checks the ordered/unordered normalization independently
at (10,4).

The exact multivariate polynomial is constructed from the continuant recurrence.
For every length 3 through 7, all variable permutations give exactly n!/2
distinct polynomials; only identity and reversal fix the polynomial. Its
homogeneous coefficient pattern, comparison degree bound and evaluations on all
words over digits 1 through 3 are checked exactly.

The source's introductory (3,1,4) illustration contains an arithmetic error:
its correct endpoints 5/19 and 6/23 give (1+5/19)/(1+6/23) = 552/551,
whereas its final logarithmic rewrite gives 551/550. Here chi(3,1,4) = 552.
The checks use the endpoints and the source's general formula (5.1), which agree;
the incorrect illustrative rewrite is excluded.

These finite checks support a proposed proof by polynomial rigidity and the
classical Schwartz–Zippel bound. They do not prove Conjecture 6, its limit or a
uniform asymptotic bound. No Lean proof is included.
