# Random-graph and nominal finite-power signatures

Source: Gabriel Ciobanu, *Finite Dependence and Invariance Hierarchies for Finitely Supported Structures*, arXiv:2609.15879v1, Question 12.3 (first target) and the proposition titled “Sort counts”, item (iii).

Run `python3 check.py` with Python 3; there are no external dependencies. The output uses exact integer arithmetic throughout and ends with `ALL_CHECKS_OK`.

The source fidelity check constructs witnesses for all fresh adjacency patterns over two context families of each size 0 through 8 in the binary-digit presentation of the Rado graph. It reproduces the source's value 2^k on these contexts. Equality-pattern enumeration and an independent Venn-multiplicity enumeration give 3 square orbits and 16 cube orbits for unordered two-element subsets. Equality and adjacency patterns give 3 square orbits and 15 cube orbits for Rado vertices. Further checks cover vertex powers through 5, support-intersection witnesses for support sizes 0 through 5, and the 7 square orbits of ordered distinct pairs.

These are finite checks of the proposed invariant. They do not prove that every nominal categorical atom with 3 square orbits is the unordered-pair object, establish the categorical invariance of the counts, or certify categorical non-equivalence. Those require the general support and category arguments in the preregistration.
