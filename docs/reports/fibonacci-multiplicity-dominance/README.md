# Fibonacci multiplicity dominance

Exact finite checks for Leonid Petrov and Jeanne Scott, *Random Fibonacci
Words via Clone Schur Functions*, arXiv:2412.21126v2, Conjecture 8.22.
Run `python3 check.py` with Python 3; no packages or input downloads are needed.
The saved result is `check.out`, ending in `ALL_CHECKS_OK`.

The program reproduces the Definition 8.10 ribbon example and the Section 8.1
partition statistics. It enumerates all set partitions through weight 9,
selects the noncrossing partitions by their arcs, and checks the source's
Catalan totals and surjectivity onto positive compositions. Every positive
composition of weights 1 through 14 satisfies the conjectured support
conditions: equal counts of twos, prefix dominance, and the displayed
lexicographic direction. All computations use integers.

Two interpretations are explicitly distinguished. `B_literal` retains
Definition 8.12 exactly, including the printed sum bound. For odd composition
length at least five its output loses the final part's weight. The program
checks the literal matrix only on its stated same-weight indices.
`B_weight_preserving` adds that final part to the first coordinate, and is a
separate interpretation. It reproduces all 64 entries of Figure 8. At weight
five the sole difference from the literal same-weight block is the entry
(221,122), which is 1 in Figure 8 and 0 in the literal block. Thus the figure's
sum is 42, while the literal block's sum is 41. This discrepancy is reported
explicitly and is not offered as a refutation of the conjecture.

Finite checks supply fidelity and regression evidence. They do not prove the
universal conjecture. The planned settlement uses run lengths and the positions
of successive twos for arbitrary compositions.
