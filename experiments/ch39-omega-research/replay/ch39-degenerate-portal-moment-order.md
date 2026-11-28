# Sharp source order and stability at a single scalar mass

Let $F$ be an $N\times N$ Hermitian portal for $N\ge1$ complex Gaussian scalars of the same known positive mass, with canonical kinetic term and known real source normalization. At any fixed set of $k\ge3$ mutually separated Euclidean source positions the propagators are scalar in the mass indices. In each oriented contraction the portal factor is $p_k=\operatorname{tr}(F^k)$. The common overall sign and normalization are known and the products of scalar propagators have positive magnitude, so division by their summed kernel recovers $p_k$. All mass-basis unitary changes are allowed. Identifying $F$ up to this equivalence amounts to identifying its real spectrum with multiplicities, including zeros. [Separated-source model](ch39-degenerate-portal-traces.md)

The shortest consecutive list beginning at order three that determines every such spectrum is $p_3,\ldots,p_{N+2}$. For spectra in $[-1,1]$, let $\varepsilon=\max_{3\le k\le N+2}|p_k(x)-p_k(y)|$ and let $\delta=\max_i |x_i-y_i|$ in increasing order. Then
$$
\delta^{N+2}\le \tfrac{2^N}{3}[4(N+1)]^{N+1}\varepsilon.
$$
The exponent $1/(N+2)$ in the inverse estimate is sharp. For an a priori interval $[-R,R]$, $R>0$, normalize the energies by $R$ and the $k$-th power sums by $R^k$. In particular, the error in the estimate is the error after the known source kernels have been divided out.

## Counting function and inverse bound

Let $h(t)=\#\{i:x_i>t\}-\#\{i:y_i>t\}$. Reading the common ordered list of endpoint positions gives $N$ upward and $N$ downward unit steps, starting and ending at zero. Ties may be expanded into separate unit steps at the same position. Every excursion of this expanded walk between successive visits to zero uses at least one step of each sign, hence there are at most $N$ excursions. The nonzero signs visible on open intervals of the original walk can therefore alternate at most $m\le N-1$ times, also for repeated or coincident endpoints. Choose one separator $r_j$ between each successive pair of opposite signs. Up to an overall sign, $Q(t)=\prod_{j=1}^m(t-r_j)$ has the sign of $h$. Integration of each indicator gives
$$
 \int_{-1}^1 h(t)t^j\,dt=\frac{p_{j+1}(x)-p_{j+1}(y)}{j+1}.
$$
The coefficient absolute sum of $Q$ is at most $2^m$. Consequently $\int |h t^2 Q|=|\int h t^2 Q|\le 2^{N-1}\varepsilon/3$.

On the open interval between a maximally displaced ordered pair $x_i,y_i$, one list has at least $i$ entries below the test point and the other at most $i-1$. Thus $|h|\ge1$ throughout this interval of length $\delta$. Remove the neighborhoods of $0,r_1,\ldots,r_m$ of radius $\eta=\delta/[4(N+1)]$. At least $\delta/2$ remains, and there $|t^2Q|\ge\eta^{m+2}\ge\eta^{N+1}$ since $\eta\le1$. Comparing the two integral bounds proves the displayed inequality, including exact injectivity at $\varepsilon=0$.

## Minimal order and sharp exponent

At a strictly ordered vector of distinct positive entries the derivative of $(p_3,\ldots,p_{N+1})$ has rank $N-1$. Its rows have entries $k x_i^{k-1}$, a Vandermonde array after removing the row factors $k$ and column factors $x_i^2$. Appending the $p_{N+2}$ row gives rank $N$. The inverse function theorem therefore gives nearby distinct positive spectra sharing all readings through $N+1$ and differing at $N+2$. For $N=1$ the shorter list is empty. Scaling both spectra by $s\downarrow0$ makes $\delta$ proportional to $s$ and the only discrepancy proportional to $s^{N+2}$. This rules out every larger uniform inverse exponent. For $N=2$, the positive spectra $(1,12)$ and $(9,10)$ have cube sum $1729$ and fourth powers summing to $20737$ and $16561$. [Integer and Jacobian checks](ch39-degenerate-portal-moment-order.py)

## Local stability at repeated eigenvalues

Fix a reference spectrum with distinct centers $a$ of multiplicities $m_a$. For two spectra in a sufficiently small neighborhood of the reference, the sharp two-spectrum inverse Hölder exponent is $1/L$, where
$$
 L=\max_a\bigl(m_a+2\mathbf1_{\{a=0\}}\bigr).
$$
Choose small disjoint cluster intervals, staying away from zero for nonzero centers. The lists have the same number of entries in each interval and their counting functions agree in the gaps. Choose internal sign separators within each cluster and any required between-cluster separators in the middle of fixed gaps. In the interval containing a maximally displaced pair, there are at most $m_a-1$ internal separators. Factors from the other separators have a fixed positive lower bound. If $a\ne0$, $|t|^2$ also has a fixed positive lower bound. The preceding deletion argument thus bounds the same integral below by a neighborhood-dependent positive multiple of $\delta^{m_a+2\mathbf1_{\{a=0\}}}$. Taking the neighborhood small enough that $\delta\le1$ gives the claim.

For sharpness in a cluster of size $m$ at $a\ne0$, choose a monic polynomial of degree $m$ with simple real roots. A small nonzero constant perturbation preserves simplicity and reality, and Newton's identities give equal first $m-1$ power sums of the two distinct root lists. Translate the roots to $a$ after scaling them by $s$. The binomial formula makes every measured difference $O(s^m)$, while the maximal displacement is proportional to $s$. For a zero cluster use the positive pair sharing $p_3,\ldots,p_{m+1}$ constructed above at size $m$. All measured differences are then $O(s^{m+2})$. Keep all other eigenvalues identical and choose a cluster attaining $L$. These pairs approach the reference and rule out any larger exponent in a neighborhood.

For comparison, Melánová, Sturmfels and Winter study power-sum recovery over several domains. Their Proposition 24 gives uniqueness for $N$ distinct positive powers on the **nonnegative** orthant, up to permutation. The signed conclusion here uses the particular consecutive powers starting at three, whose derivative exponents share the nonnegative factor $t^2$. [Recovery from Power Sums, Proposition 24](https://arxiv.org/abs/2106.13981) · [Published article](https://doi.org/10.1080/10586458.2022.2061650)
