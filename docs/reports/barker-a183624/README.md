# OEIS A183624: native-array checks

The source is R. H. Hardin's [OEIS A183624](https://oeis.org/A183624/internal),
revision 8 (22 February 2018), including Colin Barker's FORMULA conjectures.
For n >= 1, a(n) counts labelled (n+1)-by-2 arrays on {0,1,2} whose adjacent
2-by-2 blocks each sum to 4. The conjectures are

- a(n) = 2 + 2^(n+2) + 3^(n+1);
- the ordinary generating function is x*(19-69*x+54*x^2)/((1-x)*(1-2*x)*(1-3*x)).

Run `python3 docs/reports/barker-a183624/check.py` from the repository root.
The program uses Python's standard library and exact integers and fractions.
Its output is stored in [output.txt](output.txt).

The 26 exhaustive checks cover m=0,...,4 and positive row counts up to 6,
subject to an ambient-array limit of 600000. Each actual array is tested against
the original four-entry condition, with block total 2m. The largest search has
531441 arrays. All counts agree with
(m+1)^r + 2*sum(j^r for j=1,...,m).
The source's first five terms are independently reproduced by exhaustion:
19, 45, 115, 309, 859.

A second enumerator uses the actual ordered row pairs as its states and the
four-entry constraint as its transition test. For m=2 it verifies the closed
form through n=40 and the generating-function product coefficientwise through
degree 40. An exact polynomial check also verifies the numerator obtained from
2x/(1-x) + 8x/(1-2x) + 9x/(1-3x).

These finite checks provide fidelity and error detection. The unbounded
statement requires a mathematical proof via the decomposition into initial
row-sum fibers; no finite computation is claimed to prove it.
