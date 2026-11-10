# Linear coarse insolvability tests

The check concerns André Rognes, *A generalization of a representation of the
integers modulo p, for the purpose of occasionally establishing the unsolvability
of diophantine inequalities*, arXiv:2607.00757v1, Section 4.

Run with Python 3.12 or newer:

```sh
python3 docs/reports/rognes-2026-linear-coarse-strictness/check.py
```

The program uses exact integers and the Python standard library. No floating
point arithmetic, random sampling, external code or data is used. Its output is
stored in `check.out`. The final line is `ALL_CHECKS_OK`.

For fidelity it tests the sum approximation in Example 2.2 and its symmetric
intersection in Theorem 2.1 on all 16,641 pairs of integers from −64 through 64.
It checks the proposed separating input

\[
0<x<y<u<v\quad\wedge\quad x+x=v\quad\wedge\quad y+y=u.
\]

Exact addition determines u=2y and v=2x, and no input pair in that range gives
an exact solution. For moduli 2 through 17 and 40 positive multiples N of each
modulus, the tuple (N,2N,5N,6N) passes both source sum approximations and fails
the proposed doubling approximation. One common tuple also passes a finite
collection with all those moduli. The expanded approximation is checked for
soundness on every input pair and for exactness on the diagonal.

The back-and-forth checks compare all one-point extension signatures for finite
intervals: a signature consists of a residue and the two endpoint distances,
each truncated at H_j. For p=2..13, j=0..4 and interval lengths at most
H_(j+1)+4p+5, every eligible pair with equal endpoint residues has identical
signature sets, where H_0=1 and H_(j+1)=2H_j+p. It also checks 840 transfer
configurations for N=p(M+H_r+1).

These finite checks support the proposed route. They do not prove invariance for
all formulas, the unbounded rank lemma, sign elimination, or the open question.
The proposed positive settlement uses a general first-order invariance argument;
finite enumeration is not its principal mathematical content.
