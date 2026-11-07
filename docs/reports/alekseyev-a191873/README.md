# OEIS A191873: zero diagonal versus unrestricted rectangle-free matrices

Run `python3 docs/reports/alekseyev-a191873/check.py` from the repository root.
Python 3.10 or later and its standard library suffice. The stored result is
[output.txt](output.txt). The program uses integers, bit operations, and exhaustive
finite enumeration, with no floating point or external solver.

The source is Max Alekseyev's conjecture in the COMMENTS of
[OEIS A191873](https://oeis.org/A191873), revision 36, with the definitions in
[A072567](https://oeis.org/A072567) and [A001197](https://oeis.org/A001197).
A rectangle may use arbitrary distinct rows and columns. No symmetry condition
is imposed. Source definitions and their attribution remain those of OEIS;
this program independently checks the finite combinatorial facts.

The program exhausts all binary matrices for each size 1 through 4, using two
independent rectangle predicates. It reproduces the source terms: unrestricted
maxima 1, 3, 6, 9; zero-diagonal maxima 0, 2, 6, 9; and forcing thresholds 4, 7,
10 at sizes 2, 3, 4. It checks a zero permutation transversal for every enumerated
rectangle-free matrix without a full row or column, and checks every extremizer
at sizes 3 and 4. It also verifies Knuth's displayed 8-by-8 witness of weight 24
and the cyclic weight-2n family for sizes 3 through 64.

The finite checks support the source reading and proposed proof route. They do
not prove the conjecture for unbounded n, do not establish the optimality of the
8-by-8 witness, and contain no Lean proof.
