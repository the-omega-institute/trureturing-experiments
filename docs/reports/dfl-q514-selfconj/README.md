# Dietz–Flowers–Lockard Question 5.14: exact preliminary checks

Source: Geoffrey D. Dietz, Timothy B. Flowers and Shannon R. Lockard,
*A study of m-ary partitions whose conjugates are q-ary*,
arXiv:2609.08799v1; DOI: 10.13069/jacodesmath.v13i3.393,
Question 5.14, printed page 404. Source article: CC BY 4.0.
The check program is independently written; it contains no copied source code.

Run `python3 check.py`. Python 3 standard library only; all evidence uses
exact integers. The recorded run is in `output.txt`.

Fidelity checks reproduce the Definition 1.1 example and conjugate,
Proposition 1.3 on 24,992 multiplicity tuples, and the square family in
Discussion 5.11 on 189 cases. An independent cell-reflection calculation
checks every partition of each integer 1 through 35 (81,155 diagrams).

The program checks the proposed identity
`N - (2L - 1) = sum_{j<K} a_j (m^j - 1) + (a_K - 1)(L - 1)`
and divisibility of its left-hand side by `(m - 1)^2` for self-conjugate
m-ary diagrams. It also checks the odd-base consequence
`N = 1 mod 2(m - 1)` and the absence of odd power sizes in the tested family.
The compressed enumeration checks 511 nonempty choices of distinct row
lengths from `{m^0,...,m^8}` for every base 2 through 22: 10,731 diagrams,
including 5,110 diagrams at bases 3,5,...,21. Every constructed diagram is
checked against actual Ferrers transposition of its constant-width blocks.
The even-base diagram `(4,2,1,1)` is checked as a self-conjugate binary
partition of 8.

These are finite checks supporting a preregistration, not a proof of the
universal statement. No Lean proof is included.
