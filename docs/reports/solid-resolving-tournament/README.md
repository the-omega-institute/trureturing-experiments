# Solid-resolving sets in tournaments

Source: Anni Hakanen and Pavan P D, *Solid-Resolving Sets on Directed Graphs*,
arXiv:2610.35041v1, Conjecture 3.4, Definition 1.1, Theorem 2.1 and
Propositions 3.2–3.3. The conjecture asks for the six universal lower bounds
at orders 6, 7, 8, 9, 10 and 11.

Run with Python 3 and its standard library:

```sh
python3 check.py > output.txt
```

The program uses integer distances and `None` for infinity. It independently
compares Definition 1.1 against Theorem 2.1 on all tournaments and landmark
sets through order four (1,098 comparisons), reproduces Proposition 3.2 on
all labelled tournaments of orders three to five, and checks Proposition 3.3
on transitive orientations of orders two to eleven. The singleton is excluded
from the latter check: literally, Definition 1.1 makes its empty landmark set
solid because there is no nonempty vertex subset distinct from its singleton.
This boundary does not occur in Conjecture 3.4.

For the proposed general proof, it checks all 8, 64 and 1,024 labelled landmark
tournaments of orders three, four and five. In the strong case it verifies
the representative orbits, exact bounded feasible profile sets and their
explicit chain partitions. In the nonstrong case it checks the terminal
component sizes and the exact widths of the binary profile relaxations.
These are helper checks for a general antichain and terminal-component
argument, not an enumeration of all tournaments at orders seven to eleven.
A separate direct check verifies the first conjectured bound on all 32,768
labelled tournaments of order six.

Bits are assigned to pairs `(i,j)` with `i<j` in lexicographic order,
least significant bit first; a one means the arc `i -> j`. Chain data use
one-based landmark labels. Output is in `output.txt`, ending in
`ALL_CHECKS_OK`. No floating point arithmetic or external package is used.
The checks are computational evidence; no Lean proof of either direction
of the conjecture is supplied here.
