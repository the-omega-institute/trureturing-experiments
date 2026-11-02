# Individually constrained matroids and EF1

Exact computational evidence for Conjecture 4, Section 5, page 21 of
Anholcer, Bartkowiak, Bosek and Grytczuk, *Epistemic fair division of independence
structures*, [arXiv:2606.11494v1](https://arxiv.org/abs/2606.11494v1).
The target is ordinary EF1 for additive nonnegative valuations, separately for
goods and chores.

Run with Python 3, without dependencies:

```sh
python3 docs/reports/anholcer-et-al-epistemic-ef1-matroids/check.py
```

The stored output is `check.out`. All arithmetic is integral or rational.

The source fidelity check reproduces its Section 4 seven-edge graphic-matroid
example: the value matrix is [[9,2,2],[3,8,2],[3,2,2]], the displayed allocation
is EF1, every cyclic five-edge interval is a spanning tree, and the displayed
nonlinear valuation gives both bundles value 2.

For two agents, the four-item witness has capacity-one partition blocks
{1,2},{3,4} for the first matroid and {2,3},{4,1} for the second.
Both singleton-value vectors are (3,1,3,1). Both matroids are loopless, rank 2
and chromatic number 2. Exhausting all 16 labeled allocations leaves precisely
the two odd/even allocations. Their bundle values are 6 and 2, and removing
any object from the richer bundle leaves 3. Neither allocation is EF1 for goods
or chores. All 256 corners of the independent eight-coordinate radius-1/10
perturbation box also fail both notions.

The program exhausts all labeled allocations on alternating cycles of lengths
4,6,...,16. It checks rank and two-colorability, validates all independence
axioms through rank 4, tests nine strictly positive rational parameter choices
per rank, and checks equality and lower-side controls at the threshold
(r−1)a = rb. The mathematical mechanism is rank saturation: each bundle
must meet every block of both alternating matchings once, so it is a parity
class. Positive weights a on odd elements and b on even elements fail EF1
whenever (r−1)a > rb. The unbounded-family argument is a preregistration
outline; this program checks only the stated finite range.

Dror, Feldman and Segal-Halevi, *On Fair Division under Heterogeneous Matroid
Constraints*, JAIR 76 (2023), 567–611,
[DOI 10.1613/JAIR.1.13779](https://doi.org/10.1613/JAIR.1.13779),
supply the complementary-item obstruction for common nonmatroid constraints
(arXiv:2010.07280v4, Section 3.3). Their heterogeneous-category obstruction in
Section 3.2 uses a capacity-zero block, so that individual matroid cannot be
colored finitely. Their eight-item unequal-capacity illustration has individual
chromatic numbers 3 and 2. The proposed witness instead meets the individual
two-colorability hypotheses. These checks do not establish literature priority
or provide a Lean proof.
