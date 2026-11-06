# Hardin A185767: eight-column trace arrays

Run `python3 check.py` in this directory. The stored standard output is
`check-output.txt`; the program exits 0 and ends with `ALL_CHECKS_OK`.
All arithmetic used as evidence is exact Python integer arithmetic.

The objects are labelled binary arrays with eight columns and n+1 rows.
For each contiguous 2-by-2 block its trace is the sum of its main-diagonal
entries. Traces must differ at horizontally or vertically adjacent block
origins. Boundaries are open; there is no symmetry quotient.

The checker independently enumerates all 65,536 two-row arrays, validates
the three source examples, and obtains the first four counts by transposing
the arrays and transferring across a different state space. It constructs
the 4,456 valid row-pair states and 32,456 overlapping transitions directly
from the predicate. Forward path counts match all 200 source terms.

For the printed degree-19 polynomial p, direct powers and a separate Horner
evaluation agree on all eight residual vectors. The checks establish
T^7 p(T)1 = 0 in every coordinate, while the coordinate sums of T^5 p(T)1
and T^6 p(T)1 vanish. The recurrence residual at n=24 is 78,336,000; those
at n=25 through 200 vanish. The unbounded conclusion additionally requires
the array/walk bijection and propagation argument; this report supplies
arithmetic evidence for that planned proof, with no Lean settlement.

Source: R. H. Hardin, [OEIS A185767](https://oeis.org/A185767), NAME, FORMULA,
EXAMPLE and the linked [b-file](https://oeis.org/A185767/b185767.txt).
The entry is revision 9, October 3, 2025; the original submission is
February 3, 2011. The unmodified b-file is stored as `b185767.txt`.
OEIS data are subject to the [OEIS End-User License Agreement](https://oeis.org/LICENSE).
The checker is independently written for these literal definitions.
