# Toroidal grids and unweighted caterpillar PCGs

Source: C. Antony, M. Lafond, A. Monti and B. Sinaimeri, *Pairwise Compatibility Graphs of Unweighted Caterpillars*, arXiv:2610.14681v1, Problem 4, https://arxiv.org/abs/2610.14681v1.

Run with Python 3.12 or newer:

    python3 docs/reports/antony-lafond-monti-sinaimeri-2610.14681-toroidal-grids/check.py

The program uses exact integers and finite sets, with no external packages. It reproduces the paper's explicit cycle labellings for every order from 3 through 32. It constructs the corresponding unit-edge caterpillars and checks actual breadth-first path distances, all leaves, and the spine after leaf deletion.

For the proposed negative answer to Problem 4 it checks the scalar interval lemma on its stated finite range, every vertex and every pair of potential extreme neighbours of all toroidal grids with both cycle orders from 3 through 12, and all ordered products from six additional minimum-degree-two examples. The graph obstruction is checked both by searching all local witnesses and by the proposed three-case construction. Cycles provide a negative control: they admit representations and are not rejected by the obstruction.

The stored output records the exact ranges and counts. These are finite checks supporting a general proof route; they do not certify the unbounded theorem or a Lean proof. The intended settlement uses a general inequality and Cartesian adjacency argument, rather than positive exhaustive enumeration. The check is original experimental code; no third-party program or Lean settlement code was copied.
