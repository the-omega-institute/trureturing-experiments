# Independent-domination root closure: exact finite checks

Run `python3 check.py` with Python 3. The program uses only standard-library integer and rational arithmetic. The stored output is in `output.txt`; a successful run ends with `ALL_CHECKS_OK` and exits 0.

The fidelity check reproduces the displayed identity
$D_i(G[K_r],x)=D_i(G,rx)$ in Shao and Zhao, *On the Roots of Connected Domination Polynomials*, arXiv:2609.06569v5. It enumerates all labelled simple graphs on at most four vertices and all $r=1,2,3$.

The proposed settlement is supported by independent subset enumeration for every complete multipartite graph whose positive part sizes partition an integer from 1 through 12. The resulting polynomial is the sum of one monomial per part. In particular, the program checks $D_i(G(a,b,m),x)=a x^m+b x^{2m}$ for every positive $a,b,m$ with $m(a+2b)\le14$. It verifies connectivity, exact divisibility by $a+b x^m$, zero roots, the rational roots $-a/b$ when $m=1$, and the nonreal root $i$ of the polynomial of $K_{2,4}$. Separate checks enumerate all labelled graphs on at most five vertices and verify nonnegative nonzero coefficient vectors and positivity at three positive rational arguments. The empty graph has polynomial 1.

These finite checks do not prove density or either universal closure equality. The proposed universal argument uses whole-part classification, rational density and the mesh of odd root angles; it is independent of the finite bounds used here. No numerical approximation or floating-point root finding supplies evidence.

The program and this report are original contributions under the repository license. The cited source is not bundled.
