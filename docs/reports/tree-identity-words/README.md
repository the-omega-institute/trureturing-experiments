# Identity words on trees

Exact finite checks of Boris Shapiro, *Transposition Cayley Spectra and Subgraph
Counts*, arXiv:2610.14373v1 (9 October 2026), Remark 8 and Conjecture 1(iii).

Run with Python 3, without third-party dependencies:

```sh
python3 docs/reports/tree-identity-words/check.py
```

`output.txt` contains the complete output. The program exits 0 and ends with
`ALL_CHECKS_OK`. Every permutation, count and comparison uses integers.

The source's 13 tree-table values are reproduced: 1, 4, 30, 42, 336, 608, 816,
5040, 11220, 14480, 18880, 19680 and 25520, with each value matched to its named
tree. Prüfer codes and canonical rooted-tree encodings produce one representative
of every tree isomorphism class on one through six vertices: 1, 1, 1, 2, 3 and 6
shapes, respectively. Relabeling preserves both sides of the assertion.

For every representative, the program enumerates every double-occurrence word
and compares its literal permutation product with the vertexwise noncrossing
condition. There are 688,148 such words in total. An independent meet-in-the-middle
count groups all half words by permutation and used-edge mask, without assuming
double occurrence. It counts all length-2e identity words using every edge; equality
with the enumerated double-occurrence identity count rules out any additional
identity word with different multiplicities. Literal brute force over all words for
every tree with at most four edges cross-checks that counting algorithm.

Products use `(sigma*tau)(x) = sigma(tau(x))`. Distinct letters cross exactly when
their occurrence positions alternate. A crossing is forbidden precisely when the
two edges share an endpoint. Nested occurrences are allowed.

These checks reproduce the source's finite verification range. They do not prove
the statement for arbitrary trees and do not certify literature novelty.
