# Strong total cover with integer demands

Source: Manuel Aprile and Luca Ferrarini, *The Strong Total Cover Problem*,
[arXiv:2610.23312v1](https://arxiv.org/abs/2610.23312v1), Section 5's first
algorithmic question and its explicit constant-factor approximation belief.
The source suggests factor 2; these checks concern a proposed factor 3.

Run with Python 3.12 or later:

```sh
python3 check.py
```

`output.txt` contains the complete successful output. The final line is
`ALL_CHECKS_OK`. The implementation is original and uses the standard library.
All evidence is integer arithmetic or `fractions.Fraction`; no floating-point
LP solver is used.

Fidelity checks reproduce the source's unweighted complete-graph optimum for
orders 2 through 5, and its weighted integrality-gap construction for complete
graphs of orders 3 through 6. Exhaustive integral optimization and matching
rational primal/dual certificates independently check the latter LP values.

For rounding, every labelled simple graph of order at most 3 is checked with
fractional coordinates on the thirds grid in [0,1], all 0/1 cost numerators
(equivalently costs 0 and 1/2), and all feasible vertex demands. Checking the
largest feasible edge demand covers all smaller edge demands because the output
is independent of edge demands. Unrestricted multiplicities are checked on all
graphs of order at most 2, all coordinates in {0,1/3,...,8/3}, every feasible
demand vector, and every such cost vector. An additional 20,000 seeded tests
use graphs of order at most 8, costs with denominator 7, and coordinates with
denominator 12. Seed: 261023312.

The checks test edge coverage by the threshold set, the cheapest-item fractional
bound, the appearance count in vertex rows, repair feasibility, the cost bound,
and floor/residual recombination. Empty graphs, isolated vertices, zero costs,
zero demands, threshold equality, and repeated selections are included.

Finite checks do not prove the universal certificate and do not formalize
polynomial bit complexity. No claim that factor 2 is settled is made.
