# Abelian closure of finite-type shifts

`check.py` is an exact diagnostic for Karhumäki–Puzynina–Whiteland,
*On Abelian Closures of Infinite Non-binary Words*, arXiv:2012.14701v1,
Section 6 and Problem 7.3. It uses Python 3 and its standard library.

Run:

```sh
python3 -I docs/reports/karhumaki-abelian-closure-sofic/check.py
```

`output.txt` contains the output. Successful execution ends in `ALL_CHECKS_OK`.
All evidence uses integer counts, finite sets and exact reachability; there is
no floating point or randomness. These diagnostics do not certify the infinite
claim or its negation in Lean.

The graph has vertices 0 through 5 and the following non-reset edges:

| Source | Label | Target |
| --- | --- | --- |
| 4 | a | 4 |
| 4 | b | 4 |
| 5 | a | 5 |
| 5 | b | 5 |
| 0 | c | 1 |
| 0 | c | 5 |
| 1 | a | 2 |
| 2 | b | 1 |
| 1 | d | 3 |
| 4 | d | 3 |

Its 15 reset edges, labelled e, are all edges v to 0, all edges 0 to v,
and the loops at vertices 1, 2, 3 and 4. Duplicates are removed. Without e,
any path containing both c and d is exactly c(ab)^k d. The free a,b loops
at 4 and 5 supply paths with only one marker. The reset edges make the graph
strongly connected.

The intermediate alphabet consists of the 25 annotated edges and five plain
payload letters. A lift word is legal when it has no three consecutive plain
symbols and every subword of length at most four admits a graph path respecting
the annotations. The program exhausts the finite combined automaton of recent
symbols and globally reachable graph states. It finds 2,334 reachable states
and 18,725 allowed transitions, with no local/global mismatch.

For r at least 5, each lifted symbol is coded as a delimiter, its payload,
and a permutation of r distinct tags. Plain symbols use the identity tag
permutation; the annotated edges use the first 25 nonidentity permutations in
lexicographic order, in the edge order printed in the program. Codewords have
length r+2. The output alphabet has r+6 letters. Delimiter spacing, codeword
membership and the lift rules give a finite forbidden set of maximum word
length 4(r+2).

The target point has canonical blocks with payloads e infinitely in both
directions, surrounding c a^k b^k d. Any finite block hull containing both
markers can be replaced by c(ab)^k d with the same exterior e letters. A hull
with at most one marker is already a graph path. The representatives preserve
both endpoint letters. Keeping those endpoint blocks plain preserves every
partial-block trimming, while the interior annotation permutations preserve
each tag count.

The checks comprise:

- Source fidelity: exhaustive recovery of the source's ternary graph from its
  five forbidden pairs, and the golden-mean all-subfactor Parikh test on all
  8,191 binary words of lengths zero through twelve.
- Exact graph reachability and the potential certificate for the balanced
  c-to-d component.
- Exhaustion of the finite lift-memory automaton.
- All 4,260 nonempty block hulls from the displayed target windows, for
  k=0 through 8 and r=5 through 8, with all 300,510 valid endpoint cuts.
- Independent graph Parikh feasibility for all 169 pairs m,n=0 through 12,
  giving feasibility exactly when m=n, and 156 ordered distinguishing pairs.
- The two annotated loops at vertex 4 give 2^n distinct paths for n=0 through 12.

The finite cutoff checks do not establish an unbounded language identity,
soficity equivalences, the all-r statement, or an exhaustive literature screen.
Those obligations belong to the mathematical argument and subsequent proof.
