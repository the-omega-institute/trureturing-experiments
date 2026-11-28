# How hopping controls the noise amplification of two and three contact probes

Let $u>0$, $w=v/u\in(0,1)$ and $0<\rho<1$ be common to every vertex of a fixed finite simple graph with known hopping $t$. At site 1 measure $X_{AB},X_{AA},X_{BB}$, and at site 2 measure $X_{AB}$, with distinct positive degrees $z_1,z_2$. For the full parameter triple the worst local inverse amplification of errors in these four raw ratios is $\Theta(|t|^{-8})$ as $t\to0$. At a third site a single extra mixed ratio can improve the exponent to six. The criterion uses both the rooted triangle counts and the particular parameter values. Norms here are ordinary Euclidean norms on the unscaled ratios and the parameters $(\rho,u,w)$.

## The two independent low-order readings

Write $z_j$ and $n_j$ for the degree and the number of triangles through vertex $j$. If $S(\rho)$ is the three-component static curve, define $L=u^{-2}\ell(w,\rho)$ and $M=u^{-3}m(w,\rho)$. Through sixth order the three ratios at site $j$ are $S(e_j)+O(t^8)$, where
$$
 e_j=\rho+t^4z_j L+t^6 n_j M .
$$
At fourth order only a bond contributes. At sixth order a connected cluster meeting the observed site is a triangle with two hops on each of its three edges. The contribution is the same for all rooted triangles by uniform couplings. The local insertion response is parallel to $S'$, and the powers of $u$ follow from the homogeneity of the moments. The remainder holds in $C^1$ on each compact interior box by analyticity and positive static denominators. The mixed slope $S'_{AB}>0$ throughout. [Graph expansion](ch39-proton-overlap-jacobian-coordination-law.py) · [Spectral-parameter signs](ch39-proton-two-degree-signs.py)

For fixed $t$ choose a unit direction orthogonal to both $d e_1$ and $d e_2$. Every available row of the four-raw-ratio Jacobian applied to it is $O(t^8)$, so the smallest singular value is at most $C|t|^8$. The [two-degree inverse](ch39-proton-two-degree-inverse.md) gives the reverse inequality $c|t|^8$ uniformly on each prescribed compact interior box. In particular, no differentiable local inverse on the image can improve the eighth-order worst-case error exponent. The original compact-box injectivity also holds after adding an observation.

## The triangle coefficient

To extract the mixed response of a single triangle at $u=1$, let $Z_0^{ij}$ be the static one-site insertion moments, and $C^{ij}$ the three-site moments at order six with insertions only at the observed vertex. These unnormalized moments include the $3216$ allowed signed hopping configurations of the triangle. Each static three-site moment is $Z_0^{ij}(Z_0^{00})^2$. Differentiating the quotient $X_{AB}=Z^{00}Z^{11}/(Z^{10}Z^{01})$ at this static table and inserting $C^{ij}$ gives $mS'_{AB}$. The replay computes every coefficient over the rational function field in $(w,\rho)$ and checks the six intermediate moments on the previously evaluated slice $\rho=1/2$. [Exact triangle and graph replay](ch39-proton-contact-noise-orders.py)

Let $A_j=X_{AB,j}$ and choose $c=(z_2-z_3,z_3-z_1,z_1-z_2)$. Subtracting the three mixed readings with coefficients $c_j$ cancels their constant and fourth orders. If $c\cdot n\ne0$, the following normalized readings converge in $C^1$, with error $O(t^2)$:
$$
 \left(A_1,\frac{A_2-A_1}{(z_2-z_1)t^4},
       \frac{\sum_j c_jA_j}{(c\cdot n)t^6}\right)
 \longrightarrow (S_{AB},S'_{AB}L,S'_{AB}M).
$$
At fixed $\rho$, the determinant of $\partial(L,M)/\partial(u,w)$ is $u^{-6}G(w,\rho)$, where $G=3m\partial_w\ell-2\ell\partial_w m$. Wherever $G\ne0$, the limit has full rank. The raw readings then have smallest singular value at least $c_0|t|^6$ locally, uniformly on compact subsets avoiding $G=0$. A unit direction annihilating $d\rho$ and $dL$ makes every raw derivative $O(t^6)$ and proves sharpness.

If $c\cdot n=0$, the third effective shift through order six is an affine combination of the first two shifts, which leaves the two-site upper bound $O(t^8)$ on the smallest singular value. Now suppose $c\cdot n\ne0$ but $G=0$ at the parameter under study. The first two effective-shift covectors span a plane that is $O(t^2)$ from $\operatorname{span}(d\rho,dL)$: use $de_1$ and $(de_2-de_1)/[(z_2-z_1)t^4]$. These limiting covectors are independent because the two-degree inverse has full rank. At this point $dM$ lies in the limiting plane. Thus a unit vector annihilating $de_1,de_2$ gives $dM=O(t^2)$ in that direction. The residual third effective shift is $(c\cdot n)t^6 M/c_3$, yielding $O(t^8)$ on the same vector. In both exceptional cases the four original ratios still ensure the matching lower singular-value bound. These orders describe infinitesimal worst-case errors at a fixed parameter. The sixth-order test alone need not yield a global sixth-order error bound on boxes crossing its exceptional set.

## The nineteen-site shell

With adjacency defined by squared distance at most $804$, the degree, rooted triangle count and class size are $(10,45,6)$, $(14,83,3)$, $(16,100,6)$ and $(18,117,4)$. For the four increasing triples of distinct degree classes in lexicographic order the determinants of rows $(1,z_j,n_j)$ are $-8,-16,-8,0$. The triple $(14,16,18)$ stays at exponent eight. The three triples containing degree 10 attain exponent six at generic physical parameters. In fact
$$
 G(1/2,1/2)=-101786027/52502068096 .
$$
There is also an interior exception for those triples. At $\rho=1/100$ the exact rational values of $G$ have opposite signs at $w=1/100$ and $w=1/2$. Its denominators are positive between these points, so $G$ vanishes at an intermediate $w$. Two sites in the same degree class have the same rooted triangle count and give no new sixth-order contrast. [Exact forward coefficients and sign tests](ch39-proton-contact-noise-orders.py)
