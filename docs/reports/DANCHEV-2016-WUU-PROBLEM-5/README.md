# Danchev: WUU rings, Problem 5

Source: Peter V. Danchev, *Weakly UU rings*, Tsukuba Journal of Mathematics
40(1) (2016), 101–118, DOI [10.21099/tkbjm/1474747489](https://doi.org/10.21099/tkbjm/1474747489).
Problem 5 asks whether every WUU ring has abelian quotient by its Jacobson
radical. Here abelian means that all idempotents are central.

Run with Python 3, without third-party packages:

```sh
python3 check.py
```

`output.txt` contains the program output. The program uses integers modulo
positive integers, exact matrix operations, and finite power-orbit detection
for nilpotence. Inverses are constructed only after an exact coprimality test
and are verified on both sides.

The fidelity check reproduces Example 2.2(c) for upper triangular matrices of
sizes 2 and 3: over F₂ these rings are WUU; over F₃ they are not. In the WUU
examples, the program computes the Jacobson radical by testing every ring
element, recovers the strictly upper triangular ideal, and checks that every
idempotent commutes modulo that ideal.

For every idempotent and every pair of opposite-corner elements in M₂(Z/mZ),
2 ≤ m ≤ 7, and T₃(F₂), it verifies the two unit/inverse identities at t = 1, 2
and the implication from their nilpotence alternatives to nilpotence of the
corner product. It also checks the passage from are to ar by the power
identity for every listed e, a and r and exponents 0 through 3. Nonzero
nilpotent corner products occur in the checks.

These are finite diagnostics for a proposed general argument. They establish
neither the universal claim nor its literature novelty and contain no Lean
proof.
