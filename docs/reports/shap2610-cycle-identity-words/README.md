# Cycle identity words

Exact finite checks of Boris Shapiro, *Transposition Cayley Spectra and Subgraph
Counts*, [arXiv:2610.14373v1](https://arxiv.org/abs/2610.14373v1), Section 3.1 and
Conjecture 1(ii). The source reports checking the cycle formula through six
vertices. This program reproduces that range and extends it through eight.

Run with Python 3.12 or newer, without third-party packages:

```sh
python3 docs/reports/shap2610-cycle-identity-words/check.py
```

The saved output is `output.txt`. The program exits 0 and ends in
`ALL_CHECKS_OK`. All evidence uses integer arithmetic.

Each edge of the labelled cycle is a transposition; words are ordered and use
every edge. Permutation arrays contain the images of the literal product, with
the rightmost factor acting first. Appending a transposition swaps two positions
in this image array. Words are not identified under rotation or automorphisms.

The exhaustive meet-in-the-middle search inspects every half-word, groups them
by product and support mask, and joins exactly those pairs whose concatenation
has identity product and full support. It assumes neither a multiplicity bound
nor a proposed lemma. Direct Cartesian enumeration separately reproduces all
words for three, four and five vertices. For eight vertices the cycle search
examines all 2,097,152 half-words of length seven.

| Vertices | Path endpoint words | Cycle identity words | Marked singleton pairs |
| --- | ---: | ---: | ---: |
| 3 | 2 | 12 | 24 |
| 4 | 6 | 72 | 144 |
| 5 | 20 | 400 | 800 |
| 6 | 70 | 2100 | 4200 |
| 7 | 252 | 10584 | 21168 |
| 8 | 924 | 51744 | 103488 |

The path search independently enumerates every word of length `2*n-3` with
endpoint-transposition product. It checks inversion increase at every step,
one mutual extreme-token exchange, one singleton edge and all other edges
twice, and equality with the independent two-interleaving construction. Its
meeting-edge distribution equals the squares of the binomial coefficients.

For every enumerated cycle word, the program checks exactly two singleton
edges and multiplicity two for all others. Marking either singleton and rotating
it to the end produces a path endpoint word. The program checks extraction,
reconstruction, and all triples consisting of a marked edge, its original
position and a path word. It also checks the marked-pair count and the Catalan
formula independently.

These checks support the proposed general counting argument. They prove only
the stated finite ranges; they are neither a Lean proof nor an unbounded
settlement. The program is original code for these checks and uses the Python
standard library only. The external paper supplies the mathematical statement,
not program source.
