# OEIS A400412 exact arithmetic checks

Run `python3 check.py` using Python 3 with its standard library. `output.txt`
is the program output. Successful execution exits 0 with `ALL_CHECKS_OK` as
its final line. Evidence uses integers and `fractions.Fraction` throughout.

The source is [Hoft, OEIS A400412](https://oeis.org/A400412/internal),
revision 23, with its explicitly equivalent arithmetic definition: a maximal
2-dense block is a consecutive portion of the increasing positive divisor list,
cut precisely when the next divisor is greater than twice the previous one.
No restriction to odd witnesses is imposed.

The fidelity check reproduces the first nine rows from the source's search to
3025. Positive cells are exact minima; an unfilled cell in a bounded search is
never certified empty by that search.

The main checks compare every actual divisor list for two families with an
exact odd-divisor threshold formula: `2^e*q^(h-1)*r` for distinct odd primes,
`2 <= h <= 16`, `q <= 43`, `r <= 199`, `0 <= e <= 6`; and pure odd prime
powers with `2*p` odd divisors. They check the interlacing formula and the
obstruction that an even block count greater than two must be at least `h`.
These are finite checks supporting the preregistered informal proof route;
they are not a universal proof or a Lean certificate.

The family check exhausts every positive integer through 2,000,000 and records
all exact positive minima found in rows 1 through 60. Actual full divisor lists
are counted; the threshold formula is checked separately against them. Trial
division independently checks every divisor list through 20,000. The table's
`?` cells are unknown: their minimum exceeds the search bound or they are
empty. Thus the table is a partial computation of T, not a completed triangle.
The output reports any even minimum in a column at least three. Absence of
such a minimum in this finite range does not settle Conjecture (1).

Definitions and source values: OEIS Foundation Inc.,
[OEIS license and citation terms](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement).
The check program is independently written; it is not a translation of the
source's Mathematica program. No source PDF or third-party proof is copied here.
