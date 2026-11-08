# Exact checks for an all-nOA decision algorithm

`python3 check.py` writes the results in `output.txt`. The checker uses only Python's standard library, integers and Boolean tables. It exits zero with final line `ALL_CHECKS_OK`. These are finite supporting checks, not a proof of the general algorithm.

The source is M. Pavičić and N. D. Megill, *Quantum Logic and Quantum Computation*, published handbook chapter (2007), arXiv:0812.3072v1, Definition 7, equations (4)–(6), footnote 5 and the Section 8 request for an all-nOA dynamic program. The candidate applies simultaneous whole-matrix pivot updates until a complete sweep is unchanged. It tests the source inequality against that terminal matrix.

For every input, the checker constructs exact order, meet, join and complement tables and validates the ortholattice axioms. It checks inflationarity, the strict-height charge bound, terminal fixedness for every pivot, replay of one common changing-pivot word, and the terminal criterion. An independent scalar evaluator of the source recursion checks the word-to-matrix correspondence and domination on all valuations of 3OA, 4OA and 5OA until the first violating assignment; positive results exhaust all valuations. For 5OA, common scalar prefixes are cached without combining different hidden assignments. Each negative terminal answer is checked against an actual scalar-recursion violating word and a one-pivot padding.

Inputs and results:

| Input | Elements | Height in strict steps | 3OA | 4OA | 5OA | Terminal all-family criterion |
|---|---:|---:|---|---|---|---|
| Powerset Boolean algebras of sets of size 0, 1, 2, 3 | 1, 2, 4, 8 | 0, 1, 2, 3 | pass | pass | pass | pass |
| O6 | 6 | 3 | fail | fail | fail | fail |
| MO2 | 6 | 2 | pass | pass | pass | pass |
| L28, 13 atoms, 7 blocks | 28 | 3 | fail | fail | fail | fail |
| L36, 17 atoms, 10 blocks | 36 | 3 | pass | fail | fail | fail |
| L46-7, 22 atoms, 13 blocks | 46 | 3 | pass | pass | fail | fail |
| L46-9, 22 atoms, 13 blocks | 46 | 3 | pass | pass | fail | fail |

Boolean-d denotes the powerset algebra of a d-element set.

The finite fidelity results reproduce the source's O6/MO2 behavior and the authors' published OA separations. The original mathematical data for L28, L36 and both L46 examples are in arXiv:quant-ph/0009038v3 (TeX encodings next to the corresponding diagrams). The program contains those transcribed encodings. Calculations use distinct integer indices; in displayed Greechie witnesses the numeral 1 denotes the atom labelled 1, which is internally distinct from the top element. No third-party program or Lean source is used or copied.

## Source-data boundary

The later handbook chapter's literal 13-atom caption and drawing encode `123,345,567,789,9AB,BC1,BD5`. Its standard rank-three Greechie pasting is not a lattice: atoms 1 and 5 have incomparable minimal common upper bounds 3′ and B′, and no join. The checker explicitly detects this malformed input; it is outside the algorithm's ortholattice domain. The valid original L28 input is `123,345,567,789,9AB,BC1,6DC`. It is identified separately rather than silently treated as the chapter diagram.

The computed terminal results agree with the directly checked orders. The general conclusion that terminal acceptance is equivalent to every nOA law, and the universal polynomial bound, remain mathematical proof obligations. The checks do not enumerate all finite ortholattices or prove the hierarchy collapses at a uniform order.
