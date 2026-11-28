# A local conditional integral for the first Regge shape moment

Use the squared-edge Lebesgue measure, the $4^4$ periodic Freudenthal lattice and the admissible cone $C$ with simplex fatness at least $1/16$. The shape law on $D=C\cap\{V\le1\}$ is normalized by its Lebesgue volume. Put $T=S_R/\sqrt V$, $a=1920$ and let $F_\sigma=\sum_{h\subset\sigma} A_h(2\pi/n_h-\theta_{\sigma h})$, so that $S_R=\sum_\sigma F_\sigma$. The local hinge valences are four and six. [Geometry and normalization](ch39-regge-shape-support.md)

Give $C$ the probability density $e^{-V}/(\Gamma(a+1)|D|)$ in squared-edge coordinates, and call this law $\rho$. For any measurable set of shapes the volume under the cutoff $V\le v$ is $v^a$ times its volume under $V\le1$. Thus the shape marginal of $\rho$ agrees with the original one, and $V$ under $\rho$ is independently gamma distributed with shape $a$ and unit rate. For $k=1,3$ this gives
$$
 \mathbb E_D[T^k]=\frac{\Gamma(a)}{\Gamma(a+k/2)}\mathbb E_\rho[S_R^k].
$$
Translations and permutations of the four axes act transitively on the 6144 simplices. In particular $\mathbb E_\rho[V_\sigma]=1920/6144=5/16$ at any fixed simplex. To obtain the shape mean, multiply $\mathbb E_\rho[F_\sigma]$ by $6144\,\Gamma(a)/\Gamma(a+1/2)$, between $140.22610372$ and $140.22610373$. Cubing a single simplex contribution would discard the cross correlations in $S_R^3$.

## Ten variables conditional on their neighbors

Let $x$ be the ten edges of a reference simplex. Exactly 106 simplices have at least one edge in $x$. Their union contains another 200 edges $b$. For each of these 106 simplices let $w_\tau$ equal $e^{-V_\tau}$ on its positive-definite, fat admissible domain and zero elsewhere. Conditional on $b$, the mean of $F_\sigma$ is, for almost every boundary configuration,
$$
 \frac{\int F_\sigma(x)\prod_\tau w_\tau(x,b)\,d^{10}x}
 {\int\prod_\tau w_\tau(x,b)\,d^{10}x}.
$$
Indeed, after conditioning on every edge outside $x$, all other simplex factors cancel. The remaining function depends only on $b$, so conditioning just on $b$ gives the same expression. Its outer average uses the joint $b$-marginal of $\rho$.

For a one-edge update inside this integral, keep the other squared edges fixed and write $u$ for the variable edge. Its allowed set is an interval, up to endpoints. In each incident simplex the Gram matrix is affine in $u$, and the intersection of the positive-definite cone with $\det(G)^{1/4}\ge 16^{-1/4}\max_e q_e$ is convex by concavity of the determinant root. Intersect the intervals over incident simplices and require $u>0$. For almost every admissible conditioning point the interval has positive length. On it the conditional density is proportional to $\exp(-\sum_{\tau\ni u}V_\tau(u))$. Take the origin at a vertex outside the variable edge. Only the two symmetric off-diagonal entries for that edge vary. The complementary two-vector Gram block is positive definite, so Schur complementation expresses the determinant as a positive constant times $d_1d_2-(c-u/2)^2$, with $d_1,d_2,c$ fixed. The square root is concave there. The log density is consequently convex on its allowed interval. On any closed interior subinterval the exponential of the midpoint tangent to the log density and the exponential of its endpoint secant give integrable lower and upper bounds on the unnormalized conditional mass.

## Controlled cutoff for integration

Retain the box $B=\{q_e\le160\text{ for every }e\}$. Each edge occurs in at least twelve simplices, and $Q_\tau^2\le96 V_\tau$ for the maximum squared edge $Q_\tau$. Hence $B^c$ implies $V>3200$. In each simplex $A_h\le\sqrt3 Q_\tau/4$ by the triangle area bound for a given largest side. Also $|2\pi/n_h-\theta|\le2\pi/3$. Summing ten faces and applying $\sum_\sigma \sqrt{V_\sigma}\le\sqrt{6144V}$ gives $|T|<7000$.

Replace $\mathbb E_\rho[S_R^k]$ in the first formula by $\mathbb E_\rho[S_R^k\mid B]$, keeping the original gamma prefactor. The absolute error is below $10^{-124}$ at $k=1$ and $10^{-117}$ at $k=3$. To check this, let $G_s$ be unit-rate gamma of shape $s$ and put $p=\Pr(G_a>3200)$, $p_k=\Pr(G_{a+k/2}>3200)$. The discarded absolute integral after gamma normalization is at most $7000^k p_k$ by integration against $V^{k/2}$. The normalized error after conditioning is bounded by $7000^k(p_k+p)/(1-p)$. The moment generating function at $2/5$ bounds each tail by $e^{-1280}(5/3)^{a+k/2}<10^{-129}$. The box therefore controls truncation of the physical average. The sign of the remaining integral, including the exterior marginal for the mean, is open. [Exact periodic incidences and interval bounds](ch39-regge-local-moment-integral.py)
