# Exact checks for Italian cycle coefficients

For Shao and Zhao, *Counting Weight-k Italian Dominating Sets on Trees and Cycles*,
arXiv:2610.00108v1, Section 9 fourth future-work direction.
Source: <https://arxiv.org/abs/2610.00108v1>.

Run with Python 3 and its standard library:

```sh
python3 docs/reports/shao-zhao-italian-cycle-closed-form/check.py
```

The output is stored in `output.txt`. Integers and `fractions.Fraction` supply all
arithmetic used as evidence; there is no floating point calculation.

For n >= 3, d_I(C_n,k) counts functions on the labelled cycle with values 0,1,2,
weight k and neighbour-label sum at least 2 at every zero-labelled vertex. The
proposed answer concerns every 0 <= j <= 2n at weight k=2n-j. It has three disjoint
parts: zero-free deficit words, a finite sum over ordinary and exceptional gaps
between deficit zeros, and the all-exceptional periodic words. The executable
`closed_form` is the exact proposed finite sum. No rotations or reflections are
identified.

Checks and results:

- Reproduces every coefficient in the source cycle tables for n=3..12, and its
  Theorem 7.1 for j=0..5 through n=60 (347 in-domain coefficients).
- Enumerates all 88,560 ternary functions for n=3..10 using the original Italian
  inequality, checking every coefficient against the formula and a separately
  implemented boundary-conditioned weighted DP.
- Checks all 3,712 coefficients for n=3..60 against that independent DP.
- Checks the optional fixed-deficit polynomial expression at all 1,827 pairs
  with n=3..60 and 0 <= j < n.
- Checks the formula totals against the five initial values and the 53 instances
  with base n=3..55 of the repository's universal total-count recurrence:
  a(n+5)=2a(n+4)+2a(n+3)+a(n+2)-a(n+1)-a(n).
  Owner: D5/S3/Combinatorics/Graph/CycleItalianDominationRecurrence.lean in
  the-omega-institute/trureturing, source blob
  aa4c9d21c13abbb3d66b50902d3576a59f8e747e.

These finite checks support the preregistered universal combinatorial argument;
they do not prove it or certify literature novelty. No Lean proof of the formula
or its negation is included. The planned positive settlement uses a general
marked-zero gap bijection, weak compositions, slot choices and double counting.

The check program is independently written; numerical table values are quoted
from the cited source as factual comparison data. No source implementation or
paper prose is copied into the program.
