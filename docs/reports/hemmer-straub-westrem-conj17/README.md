# Hemmer–Straub–Westrem Conjecture 17: exact finite checks

Reference: David J. Hemmer, Armin Straub and Karlee J. Westrem, *New Identities in
the Character Table of Symmetric Groups involving Riordan Numbers*,
arXiv:2509.02796v1 (2 September 2025), Conjecture 6.1, equation (6.1);
Electronic Journal of Combinatorics 33(1), P1.57 (2026), Conjecture 17,
DOI [10.37236/14401](https://doi.org/10.37236/14401).

Run with Python 3, using only its standard library:

```sh
python3 -I check.py
```

`output.txt` contains the complete output. The program exits 0 and ends with
`ALL_CHECKS_OK` when every check succeeds. All arithmetic is integer arithmetic
or `fractions.Fraction`; no floating point is used.

The fidelity checks reproduce the source's Ev((3,2,2)) multiset, its R3(10)
example, all twenty character entries in Section 2.1, Table 1, and both totals
48 for lambda=(1,1,1,1). Ev multiplicities are retained.

The program enumerates all 272 partitions of sizes 0 through 12, computes both
weighted character sums directly, and checks the asserted even/odd formula.
The character calculation extracts the three-variable Frobenius coefficient;
every rectangular character is independently checked using two variables.
For sizes 0 through 7, 1003 distinct character pairs are also checked by a
separate Murnaghan–Nakayama rim-hook recursion. The sum of the weights 1/z_lambda
is checked to equal 1 in every tested size.

The resulting averages for n=0,...,12 are
`1,0,3,0,6,0,10,0,15,0,21,0,28`.
The n=0 value is a generating-series convention outside the published positive
integer quantifier. These finite checks supply source fidelity and sanity
evidence; they establish no universal theorem and contain no Lean proof.
