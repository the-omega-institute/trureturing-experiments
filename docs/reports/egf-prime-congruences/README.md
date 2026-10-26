# Tebni A350759: both prime congruences

Run `python3 check.py`. Python 3's standard library suffices. The committed
`output.txt` is the program's complete output; success ends in `ALL_CHECKS_OK`.
All mathematical arithmetic uses integers or `fractions.Fraction`.

The source is Mélika Tebni, OEIS A350759, revision 17 (12 March 2026),
https://oeis.org/A350759 and
https://oeis.org/search?q=id:A350759&fmt=text.
The source attributes the sequence to Tebni, 14 January 2022.
Its NAME is the signed first-kind Stirling transform of A345652, and its EGF is
`exp(-1 + (1+x)*(1-log(1+x)))`. The two COMMENTS conjectures ask, for every prime p,
that p divide a(p)+a(p+1)-2 and that p(p+1) divide a(p+1)-1.

The program reproduces all 24 source DATA values and the example a(9)=1000.
It independently compares the source integer recurrence with the exact formal
exponential through degree 220, with the signed Stirling-transform NAME through
degree 60, and with the rising-factorial expansion through degree 60. It checks
the proposed factorial-factorization identity and divisibility of every
nonconstant summand in that finite expansion, the proposed all-index congruence
through n=220, and both target congruences for every prime p<=211.
These finite checks support fidelity and reject small failures; they are not
proofs of the universal assertions. No Lean proof is included.

The rising-factorial expansion for A176118 is credited to Paul D. Hanna in that
entry's FORMULA (13 November 2023):
https://oeis.org/search?q=id:A176118&fmt=text.
The parent recurrence for the Stirling-transform check is from
https://oeis.org/search?q=id:A345652&fmt=text.

The quoted sequence data, definitions and recurrence are attributed to the
named contributors and The OEIS Foundation Inc.; the OEIS End-User License
Agreement applies: https://oeis.org/LICENSE. The checking implementation is
independently written for this experiment; no external source code is copied.
