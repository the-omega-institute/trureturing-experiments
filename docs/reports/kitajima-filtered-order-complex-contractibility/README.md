# Filtered order-complex contractibility

Source: Yoh Kitajima, *Filtered order complexes and magnitude homology of finite graded posets*, arXiv:2606.15241v2, Definition 2.11 and Remark 2.13.

Run `python3 check.py` with Python 3; no third-party package is required. The stored output is `output.txt`. All operations use integers, finite sets and bit vectors over F₂.

The fidelity check reproduces the chain example (Example 2.12) for ranks 1 through 6 and every filtration index from 1 through the rank. Each complex admits an explicitly verified elementary collapse to one vertex. The index-zero stage is discrete and is excluded from this contractibility statement, as in the source's later explicit range for chains.

For each t from 2 through 6, the program constructs K = ({0,1},{1,2}) * D_t, its subcomplex B = {0,1} * D_t, the nonempty face poset C, its upper filter U = C minus F({0,1}), and its lower ideal V = F(B). It adjoins minima a,x and maxima b,y, with a below C below b, x below exactly U, and y above exactly V, taking transitive closure. The resulting pure rank-four poset has 9+6t elements.

The checks verify transitivity, strict rank increase, all cover ranks, extrema, the comparable order-map contraction certificates for C,U,V, the t components of U intersect V, every rank-window identity, and all cone-union intersections used in the mathematical argument. They enumerate every chain, verify that boundary squared is zero, and compute exact F₂ boundary ranks. The index-three complex is independently collapsed to a vertex; the index-four complex has second Betti number t−1 over F₂. Nonzero homology over this field obstructs contractibility.

The program is a finite cross-check and supplies no Lean kernel evidence. Its five tested parameters do not prove the unbounded family statement; that statement requires the order-homotopy and cone-gluing argument in the preregistration. Acyclicity alone is never used as a certificate of contractibility.
