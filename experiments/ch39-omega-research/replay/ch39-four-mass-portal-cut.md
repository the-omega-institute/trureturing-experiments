# A four-layer Hermitian portal invisible to static spectra and two- and three-source responses

Consider four complex scalars with identical quantum numbers, a canonical Euclidean kinetic term and the one-loop effective action proportional to $\operatorname{Tr}\log(-\partial^2+D+s(x)F)$. Use the same local subtraction in both models. In mass units with $D=\operatorname{diag}(1,2,3,4)$ the two portals are

$$
F_\sigma=\begin{pmatrix}
5&\sigma i&1&1\\
-\sigma i&6&1&1\\
1&1&7&i\\
1&1&-i&8
\end{pmatrix},\qquad \sigma\in\{+1,-1\}.
$$

Every off-diagonal coupling is nonzero. The Hermitian matrices are positive definite by strict diagonal dominance. Since the squared masses are distinct, unitary changes of the mass basis that fix $D$ are diagonal rephasings. The real cycle product $F_{12}F_{23}F_{34}F_{41}=-\sigma$ survives both rephasing and simultaneous complex conjugation. The two portals are inequivalent.

## The constant pencil and lower-order insertions

Put $x_j=\lambda-D_{jj}-sF_{jj}$. The cross block between the first two and the last two indices in $\lambda I-D-sF$ equals $-s\mathbf{1}\mathbf{1}^{T}$, where $\mathbf{1}=(1,1)^T$. For diagonal blocks $A,B$, the Schur complement and the rank-one determinant formula give
$$
\det A\det B-s^2(\mathbf{1}^{T}\operatorname{adj}(A)\mathbf{1})
                    (\mathbf{1}^{T}\operatorname{adj}(B)\mathbf{1})
=(x_1x_2-s^2)(x_3x_4-s^2)-s^2(x_1+x_2)(x_3+x_4).
$$
The identity extends across singular diagonal blocks by polynomial equality. The entire characteristic pencil is independent of $\sigma$, including with independent real diagonal mass shifts. Thus every constant-background spectral observable agrees wherever it is defined, including a fixed-subtraction one-loop potential.

At zero background the propagator is diagonal in mass labels. The two-source trace at arbitrary Euclidean momentum is a sum of $|F_{ij}|^2$ times the scalar bubbles for squared masses $D_{ii},D_{jj}$. The coefficients coincide term by term. At three insertions reversal pairs the two orientations of each trace with a real kernel, leaving $\operatorname{Re}(F_{ij}F_{jk}F_{ki})$. The repeated-index terms depend only on diagonal entries and edge moduli. Every triangle on three distinct vertices uses precisely one imaginary edge and hence has zero real part. All CP-even three-insertion traces with real mass-dependent weights therefore agree, in particular the zero-background scalar triangle for every external momentum.

At fourth order the mass-weighted invariant $\operatorname{Tr}([D,F]^4)$ distinguishes the two signs, taking $848$ for $+$ and $736$ for $-$. In the mass basis each factor contributes $(D_{ii}-D_{jj})F_{ij}$, so the statistic weights quadrilaterals by their mass differences. An actual fourth derivative of the one-loop action instead weights the edges by propagators.

## A separated four-source response of the continuum probe

In four Euclidean dimensions the diagonal propagator between distinct points at distance $r$ is
$$
g_j(r)=\int_0^\infty \frac{e^{-D_{jj}u-r^2/(4u)}}{(4\pi u)^2}\,du
      =\frac{\mu_j K_1(\mu_jr)}{4\pi^2r},\qquad \mu_j=\sqrt{D_{jj}}.
$$
At four distinct source points the fourth derivative at $s=0$ follows from the $-\operatorname{Tr}(GsF)^4/4$ term of the logarithm. Cyclic rotation identifies four permutations of the sources, leaving six oriented loops, or three loops paired with their reversals. For a loop with successive lengths $r_1,\ldots,r_4$ each assignment of mass labels contributes the common overall factor times
$$
2\operatorname{Re}(F_{ij}F_{jk}F_{k\ell}F_{\ell i})
 g_i(r_1)g_j(r_2)g_k(r_3)g_\ell(r_4).
$$
Coincident labels give products of edge moduli or real diagonal entries times triangle products. Their real parts agree. Contact counterterms vanish at four distinct points.

Place the sources on a line at $L(0,1,3,7)$, $L>0$. For the remaining assignments the four masses are all distinct. For any fixed loop, $g_j(Ld)$ is asymptotic to a positive multiple of $L^{-3/2}\exp(-L\mu_j d)$. Strict rearrangement uniquely minimizes the decay exponent by assigning the largest mass to the shortest edge and continuing in this order. The three possible loops have the following minimum exponents.

| Cyclic order in units of $L$ | Sorted edge lengths in units of $L$ | Minimum exponent divided by $L$ |
| --- | --- | --- |
| $(0,1,3,7)$ | $(1,2,4,7)$ | $9+2\sqrt3+4\sqrt2$ |
| $(0,1,7,3)$ | $(1,3,4,6)$ | $8+3\sqrt3+4\sqrt2$ |
| $(0,3,1,7)$ | $(2,3,6,7)$ | $11+3\sqrt3+6\sqrt2$ |

The first minimum is strictly the smallest. In that cyclic order its edge indices are $(4,3,2,1)$, with cycle product $-\sigma$. Its reversed orientation adds the same real weight. All other contributions to the difference decay strictly faster, so the connected four-source responses differ for all sufficiently large $L$. This uses access to four separately located insertions in the canonical continuum probe. The static spectrum, bubble and scalar triangle already specified sample lower-order information. The full set of separated three- and four-source readings also identifies general four-layer portals, as shown below.

[Exact matrix and loop checks](ch39-four-mass-portal-cut.py) · [Large-argument expansion of the modified Bessel function, DLMF 10.40.2](https://dlmf.nist.gov/10.40.E2)

## Recovering an arbitrary four-layer portal from separated insertions

Let four known positive masses be pairwise distinct. For canonical Euclidean kinetic terms in four dimensions, the zero-background connected three- and four-source functions at separated points determine every Hermitian portal up to diagonal unitary rephasing and simultaneous complex conjugation. Finite samples suffice. The coefficients and normalization of the source coupling and the one-loop prefactor are taken as known.

**Separate the loop order.** Take a regular tetrahedron of edge length one. Each of its six squared distances can vary independently in $[1,1.04]$. With the first vertex as origin, its three-dimensional Gram matrix has diagonal entries equal to the corresponding squared lengths and off-diagonals $(d_{0i}+d_{0j}-d_{ij})/2$. At unit edges its least eigenvalue is $1/2$. On the stated box the maximum absolute row sum of the perturbation is at most $0.04+2(0.04)=0.12$, so the perturbed Gram matrix remains positive definite. Any face gives the corresponding triangle construction.

Choose the four successive edges of a single cyclic order. Hold the two unused squared lengths at one. In each used variable form an adjacent difference between successive values in the list $1,1.01,1.02,1.03,1.04$. Taking all four differences kills each of the other two loop orders because each omits a varied edge. Divide by the known signed loop factor. The resulting four-index array is the product of four changes of basis, one in each index, applied to $\operatorname{Re}(F_{ij}F_{jk}F_{k\ell}F_{\ell i})$. The change-of-basis matrix has entries $g_i(\sqrt{x_b})-g_i(\sqrt{x_{b-1}})$. For the triangle there is one unoriented order. Its raw three-index array on four squared-length choices is the tensor transform of $\operatorname{Re}(F_{ij}F_{jk}F_{ki})$ with entries $g_i(\sqrt{x_b})$.

**Invert the radial samples.** Write $G_i(x)=g_i(\sqrt{x})$ and order the masses by increasing square. Both $G_i$ and $-G'_i$ have positive integral weights in
$$
\int_0^\infty e^{-D_{ii}u}e^{-x/(4u)} a(u)\,du,
\qquad a(u)=(4\pi u)^{-2}\quad\hbox{or}\quad (4u)^{-1}(4\pi u)^{-2}.
$$
For ordered distinct $x_b>0$ the determinant of this matrix of integrals is the integral over $u_1<\cdots<u_4$ of the product of the two exponential determinants times $\prod a(u_j)$. This follows by expanding the determinant and permuting the integration variables. Both exponential determinants have constant nonzero sign on this chamber. Indeed a nonzero combination of $n$ exponentials with distinct real exponents has at most $n-1$ distinct real zeros. Divide by the first exponential, differentiate and use Rolle's theorem inductively. Thus their evaluation determinant never vanishes and its sign is constant on the connected ordered domain, also when the exponents vary in strict order. Absolute convergence follows from $x_b,D_{ii}>0$. The same argument for $-G'_i$ gives constant nonzero sign at any ordered four points. Integrating its rows over four successive intervals preserves this determinant sign by multilinearity. This proves invertibility of the difference matrix and of the raw triangle matrix for the fixed grid, including masses with coincident pair sums. At most $5^4$ four-source and $4^3$ three-source configurations suffice.

**Recover the phases.** A repeated triangle with indices $iii$ gives the real cube $F_{ii}^3$. The quadrangle $ijij$ gives $|F_{ij}|^4$, determining the support and every edge modulus. Rephase a spanning tree of each connected component so that its edges are positive. If the support has at most one independent cycle in total, its length is three or four and its measured real product fixes its phase up to conjugation. This includes a triangle with a leaf or an isolated vertex and a chordless quadrangle. A forest has no remaining phase.

With four vertices a support containing at least two independent cycles has five or six edges and a vertex connected to the three others. Use those three incident edges as the tree. For each remaining edge $ij$, the triangle through the central vertex determines $\cos\phi_{ij}$ of its oriented phase. Any two remaining edges meet at a vertex. Orient them in succession and close the path through the central vertex. The normalized real quadrangle product is $\cos(\phi_{ij}+\phi_{jk})$. Subtracting it from $\cos\phi_{ij}\cos\phi_{jk}$ gives $\sin\phi_{ij}\sin\phi_{jk}$, with the known sign for any reversal of the chosen orientation. Together with the squares these values determine the real vector of sines up to one overall sign, including zero components. That sign is precisely simultaneous complex conjugation. Reversal of each loop and cancellation of tree phases ensure that both allowed transformations preserve every measurement. The counterexample above establishes that the fourth source order is necessary in the worst case. This gives exact-data identification. Quantitative noise bounds depend on mass spacing and phase degeneracies.

For masses $1,2,3,4$ the determinants at the stated grid are enclosed near $5.486\times10^{-23}$ for the first four raw rows and $3.963\times10^{-30}$ for the four adjacent differences. The small determinants illustrate why exact-data inversion requires a separate noise analysis. The finite phase check uses quarter-turn phases, all 64 support graphs, and 191 rephasing and conjugation classes. [Radial and phase checks](ch39-four-source-portal-inverse.py)
