# Tree-shift spectrum: exact checks

Source: Jeremiah Buckley and Daniel Seco, *Rota universality of shifts on trees*, arXiv:2610.40089v1, Theorem 1 and the annular question in Remarks (iv).

Run `python3 check.py`. The program uses only integers and rational arithmetic from the Python standard library. Successful output ends with `ALL_CHECKS_OK`; the stored output is in `output.txt`.

The fidelity check enumerates every recursively labelled rooted tree with at most seven vertices (874 trees). With the source's backward shift matrix S and its transpose T, it verifies that S^m T^m is diagonal with entries equal to the m-th descendant counts. Thus its top eigenvalue equals the largest such count, reproducing the finite version of the source's norm-square identity.

Independent Gauss–Jordan inversion of T-zI at z=5/4 and z=-3/2 verifies all 11,826 inverse columns against the proposed descendant coefficient formula, zero ancestor coefficients and the resulting generation norm bound. Finite trees have leaves and their matrix spectra are not spectra of the infinite admissible trees.

The infinite example has one vertex o with a unary ray and a full binary subtree as its two child branches. The rootless version adds ancestors a_j, with parent(a_1)=a_2 and child(a_1)=o. Generation counts are 1 on the unary branch and 2^m on the binary branch; at o they are 1+2^(m-1) for m>0. At a_j they are 1 for m≤j, and 1+2^(m-j-1) otherwise. These formulas give the growth radii R_minus=1 and R_plus=sqrt(2). Prefix-free binary words 1^n0 identify infinitely many branching vertices with disjoint descendant sets.

The program enumerates generations through depth eight for ten seed vertices, checks the exact extremal profiles through depth 80, and solves a finite window with six unary ancestors. At z=5/4 the proposed binary-generation norm-square lower bounds increase by the exact ratio 32/25 and exceed 10^8 at depth 80.

These checks test the formulas and truncations. They neither prove the infinite spectrum assertion nor certify the infinite-ancestor vanishing argument. No floating-point approximation or finite-truncation eigenvalue is used as evidence for the infinite spectrum.
