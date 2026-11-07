# Roman equality under degree-class complementation

Exact finite diagnostics for Pinhe Chen, *Refined Nordhaus–Gaddum-Type Bounds
for Roman and Total Domination on δ-Complement Graphs*, Mathematics 14(11),
1797, DOI [10.3390/math14111797](https://doi.org/10.3390/math14111797),
Theorem 3 and the Section 8 open question on its equality cases.

Run with Python 3.10 or newer; no third-party packages or network are needed:

```sh
python3 check.py
```

The program resolves its data paths relative to itself. `output.txt` is its
complete output. It uses integers and `fractions.Fraction` exclusively.
For each set S of vertices labelled 2, the cheapest Roman completion labels
exactly V minus N[S] by 1. Minimizing 2|S| + |V minus N[S]| enumerates every
possible choice of S. An independent ternary-label enumeration checks this
minimum on every labelled graph of orders 1–5 and on the indicated small
unlabelled graphs and δ-complements.

The checks reproduce Chen's Example 2 (both Roman numbers equal 7), the
26 equality graphs in Table 7, its counts by order (1,3,3,6,5,8), and the
Table 6 mean refined slack (exactly 21498/2719, rounding to 7.91). They test
both directions of the proposed structural classification on all 13,595
unlabelled graphs of orders 3–8: equality holds exactly when every edge joins
equal original degrees and every original-degree class of size at least 3
induces C5, a disjoint union of edges, or the complement of such a union.
A separate test by connected components checks the explicit disjoint-union
family. All 1,099 labelled graphs of orders 1–5 also satisfy the classification.
There are 9,261 checks of the ordinary Roman predicate against the existing
D5 owner's k=1 active-neighbour predicate on ternary labels.

`graph3.g6` through `graph8.g6` are the unlabelled simple graph datasets from
[Brendan D. McKay's graph collection](https://users.cecs.anu.edu.au/~bdm/data/graphs.html),
retrieved from `https://users.cecs.anu.edu.au/~bdm/data/graphN.g6`, with N=3,…,8.
The program checks their SHA-256 hashes and counts before decoding. The
unlabelled coverage relies on that collection's enumeration; the labelled
orders 1–5 are generated independently by all edge subsets. Graph6 encoding
is documented on [the collection's formats page](https://users.cecs.anu.edu.au/~bdm/data/formats.html).
The code was written independently and does not copy Chen's supplementary code.

These computations are finite diagnostics. They do not prove the universal
classification, the D5 bridge, or any statement in Lean.
