# Sharp rank thresholds for cyclic Fourier probabilities

Let a density matrix on the cyclic group of order $d$ have coordinate probabilities $p$ supported on $S$. Assume distinct ordered pairs of different points of $S$ have distinct differences modulo $d$. Write $q$ for the probabilities of the Fourier basis $f_j(x)=e^{2\pi ijx/d}/\sqrt d$, and $u_j=1/d$.

## Minimum at fixed coordinate probabilities

Positivity forces $\rho_{xy}=0$ whenever either coordinate probability vanishes. Fourier inversion gives $\sum_j q_j e^{2\pi ija/d}=\sum_{x-y=a}\rho_{xy}$. Each nonzero coefficient therefore reads at most one matrix entry. Parseval gives

$$
d\|q-u\|_2^2=\operatorname{tr}(\rho^2)-\sum_i p_i^2.
$$

Sort the positive probabilities as $p_1\ge\cdots\ge p_k$. For $1\le r<k$, choose $c$ so that $\sum_{i=1}^r\max(p_i,c)=1$. The radius $R_r(p)$ from $u$ to the nearest readout of rank at most $r$ satisfies

$$
d R_r(p)^2 = \sum_{i=1}^r\max(p_i,c)^2-\sum_{i=1}^k p_i^2.
$$

Indeed, the decreasing spectrum $\lambda_*$ with entries $\max(p_i,c)$ up to $r$ and zeros thereafter majorizes $p$. The finite Schur–Horn theorem gives a positive matrix with that spectrum and diagonal. For any feasible spectrum, write $L_j$ for its first $j$ entries summed and $m=\#\{i:p_i>c\}$. Schur–Horn bounds $L_j$ below by $\sum_{i\le j}p_i$ up through $m$. Its prefix sums are concave and reach $1$ at $r$, so from $m$ through $r$ they also dominate the straight segment of $\lambda_*$. Every feasible spectrum thus majorizes $\lambda_*$ and has at least its squared norm. For a uniform diagonal this gives $R_r=\sqrt{(1/r-1/k)/d}$. With at most $r$ positive sites, the coordinate diagonal matrix attains $R_r=0$. The minimum-purity spectrum also occurs in the prescribed-norm frame potential problem. [Casazza et al., *A Physical Interpretation of Tight Frames*](https://doi.org/10.1007/0-8176-4504-7_4)

If $\|y-q\|_2\le\epsilon$ and $\|y-u\|_2+\epsilon<R_r$, every compatible state has rank at least $r+1$. At equality a minimizing state can yield an observation on the segment between its readout and $u$. For a specified probability vector, physical compatibility also requires zero Fourier coefficients at unused differences and positive semidefiniteness of the matrix reconstructed from its remaining coefficients and $p$.

## Uncertain diagonal and support leakage

Consider first the closed simplex on $S$, allowing zero mass at individual sites. If the true diagonal lies within total variation $\delta$ of $p$, the exact minimum radius is $R_r(p^\star)$. To construct $p^\star$, transfer up to $\delta$ of the smallest masses to the largest entry until the budget is spent or the largest reaches one. Its $j$-th prefix is $\min(1,\delta+\sum_{i\le j}p_i)$ for $j\ge1$. Every probability vector in the ball satisfies this prefix bound after sorting, by the total variation bound on subsets.

To see why concentrating the diagonal minimizes the radius, write $g_a(s)=s^2-(s-a)_+^2$. For any feasible spectrum $\lambda$ padded to $r$ entries, expanding its squared norm about $a$ and applying majorization to $(s-a)_+^2$ gives

$$
\|\lambda\|_2^2-\|p\|_2^2\ge 2a-ra^2-\sum_i g_a(p_i).
$$

Equality holds for $(\lambda_*,c)$, or at $a=0$ for a diagonal with at most $r$ positive sites. Thus $dR_r^2$ is the maximum of the right side over $0\le a\le1$. Each $g_a$ is convex, so this maximum decreases under concentration. The extremal diagonal belongs to the ball and the associated spectrum is attainable.

For leakage $t<1$ outside $S$, let $p$ now denote the conditional diagonal on $S$. The normalized compression has rank at most that of the original state. Projecting a purification onto $S$ and normalizing gives overlap $\sqrt{1-t}$ with the original purification. Trace distance contracts under partial trace and measurement, so the two Fourier distributions differ by at most $\sqrt t$ in total variation and at most $\sqrt{2t}$ in Euclidean norm. Thus every rank-at-most-$r$ state satisfies $\|q-u\|_2\ge R_r(p)-\sqrt{2t}$. In dimension three the pure state $(2,2,-1)/3$ has leakage $1/9$ from $\{0,1\}$, conditional probabilities $(1/2,1/2)$, and vanishing nonzero cyclic correlations, hence uniform Fourier probabilities.

## Exact examples

| $d$ | $S$ | $p$ in support order | $dR_r^2$ in ascending order of $r<k$ |
|---|---|---|---|
| 7 | $0,1,3$ | $1/2,1/3,1/6$ | $11/18,\ 1/9$ |
| 13 | $0,1,3,9$ | $1/2,1/4,1/8,1/8$ | $21/32,\ 5/32,\ 1/32$ |

The certificate constructs a positive witness for each entry and checks its diagonal, rank and Fourier deviation. In dimension four the support $\{0,2\}$ has a repeated nonzero ordered difference. The pure state $(1,0,i,0)/\sqrt2$ gives uniform Fourier probabilities. Run the exact arithmetic checks with Python and SymPy. [Replay script](ch08-fourier-rank-radius.py)
