# Sharp insertion orders for two scalar masses with a small multiplicity

Use the separated-source scalar portal with two known distinct positive mass values and known multiplicities summing to $N$. Write $P$ for the mass projection of smaller rank, $Q=I-P$, and $F$ for the Hermitian portal. At order $k\ge3$ the full position dependence yields the real traces of all cyclic words with $k$ occurrences of $F$, separated by either $P$ or $Q$. Summing either label inserts $I$. The order below counts consecutive source orders beginning at three. [Position-space extraction and equivalence](ch39-degenerate-portal-traces.md)

## One-dimensional mass block

Write the portal in blocks as $F=\left(\begin{smallmatrix}a&h^*\\h&C\end{smallmatrix}\right)$, with $m=N-1$. The all-$P$ cube gives $a^3$, fixing the real scalar $a$. The all-$Q$ words give $\operatorname{tr}(C^k)$ and determine the spectrum of $C$ by order $m+2$. A word that visits $P$ once and has $j+1$ consecutive visits to $Q$ gives $h^*C^jh$ at order $j+2$, for $1\le j\le m-1$. The alternating four-letter word gives $\|h\|^4$, hence also $h^*h$. For the distinct eigenvalues $\lambda$ of $C$, the first as many moments as eigenvalues determine $\|E_\lambda h\|^2$ by Vandermonde inversion. Within each eigenspace all vectors of that length are unitarily equivalent. This includes a vanishing vector and arbitrary multiplicities.

Thus orders through $\max(4,N+1)$ suffice. For $N\ge3$, take $h=0$ and keep $a$ fixed. The [single-mass order lower bound](ch39-degenerate-portal-moment-order.md) for the complementary block forces $N+1$. At $N=2$ use $a=C=0$. All cubic words vanish, but the alternating quartic distinguishes different $|h|$. These bounds are sharp under mass-preserving unitary equivalence.

## Two-dimensional mass block

For rank $P=2$ and $N\ge4$, the sharp worst-case order is $K=\max(6,2N-4)$. The data identify $F$ up to a unitary commuting with the known mass matrix and a possible global complex conjugation. Write $F=\left(\begin{smallmatrix}A&B\\ B^*&C\end{smallmatrix}\right)$. The pure-block words recover the spectra of $A$ and $C$ using the single-mass bound, since their dimensions are $2$ and $N-2$ and $K\ge N$. The alternating words of orders four and six give $\operatorname{tr}(BB^*)^2$ and $\operatorname{tr}(BB^*)^3$, where the powers are taken before the trace. These determine the two nonnegative eigenvalues of $BB^*$. At fixed positive sum of their squares, parametrize the ordered pair as $r(\cos t,\sin t)$ with $0\le t\le\pi/4$. Its sum of cubes decreases strictly on this interval. If the square sum vanishes both eigenvalues are zero.

In particular $\operatorname{tr}BB^*$ is known. The complete one- and two-insertion words are now available directly from the traces of $A,C,A^2,C^2,BB^*$. Their sums give $\operatorname{tr}F$ and $\operatorname{tr}F^2$. Summing measured words over the mass labels gives $\operatorname{tr}F^k$ for $3\le k\le N$, so Newton's identities recover the full spectrum with multiplicities.

Let $M_j=PF^jP|_{\operatorname{ran}P}$, $M_0=I_2$, and denote the number of distinct eigenvalues of $F$, now including zero when present, by $s$. The first $s$ matrices $M_0,\ldots,M_{s-1}$ determine all subsequent ones by the recurrence from the square-free spectral polynomial. Their traces and the inner products $G_{jk}=\operatorname{tr}(M_jM_k)$ are real word traces at degrees $j$ and $j+k$. Within the first $s$ indices the only pairs that can exceed $K$ are $(s-2,s-1)$ and $(s-1,s-1)$, and only when $s=N$. For $R(z)=P(zI-F)^{-1}P|_{\operatorname{ran}P}$, complementary minors and the two-dimensional trace identity give
$$
 \operatorname{tr}(R(z)^2)=(\operatorname{tr}R(z))^2
      -2\frac{\det(zI-C)}{\det(zI-F)} .
$$
The right side is known, since $\operatorname{tr}R(z)=\sum_{j\ge0}\operatorname{tr}(M_j)z^{-j-1}$ and its coefficients follow from the recurrence. The coefficient at $z^{-t-2}$ on the left is $\sum_{j+k=t}G_{jk}$. For each missing total degree, in increasing order, first reduce any index $\ge s$ by the recurrence, which strictly lowers its total degree. The first anti-diagonal has the missing pair twice and the second contains one missing square. Subtraction gives both entries.

Subtracting the scalar parts of the $M_j$ gives traceless Hermitian $2\times2$ matrices, or real Pauli vectors. Their Gram matrix fixes them up to a common orthogonal transformation. Rotations act by unitary conjugation and complex conjugation supplies a reflection. Interpolation now recovers the compressed spectral projections $A_\lambda=PE_\lambda P|_{\operatorname{ran}P}$ up to that same action. To see why this identifies the original portal, the map from $\operatorname{ran}P$ into the $\lambda$-eigenspace by $E_\lambda$ has Gram matrix $A_\lambda$. Equal Gram matrices give an isometry on its image, extendible to a unitary between the eigenspaces of the known equal multiplicities. Summing these unitaries matches $F$ and carries $\operatorname{ran}P$ to itself.

## A pair whose first difference is at the cutoff

In dimensions four and five, put $A=C=0$ and choose the two singular-value pairs of $B$ as $(1,0)$ and $(2^{-1/4},2^{-1/4})$. Only alternating words at even orders can survive. Both fourth-order traces equal one, and the sixth-order traces are $1$ and $1/\sqrt2$. Thus order six is necessary in both dimensions.

For $N\ge6$ let $g(z)=\prod_{a=1}^N(z-a)$ and let $Z=\operatorname{diag}(1,-1)$, $X=\left(\begin{smallmatrix}0&1\\1&0\end{smallmatrix}\right)$. Set $B_0(z)=(z^2+1)Z$ and $B_1(z)=(z^2-1)Z+2zX$. For $b=0,1$ form the real matrices
$$
 A_a^{(b)}=c\left(\frac{a^2+1}{|g'(a)|}I+\frac{B_b(a)}{g'(a)}\right),
 \qquad c^{-1}=\sum_a\frac{a^2+1}{|g'(a)|}.
$$
Each $B_b(a)$ has eigenvalues $\pm(a^2+1)$, so each $A_a^{(b)}$ is positive semidefinite of rank one. The partial fractions of $B_b(z)/g(z)$ have residues $B_b(a)/g'(a)$ whose sum is zero. Thus the $A_a^{(b)}$ sum to $I_2$. Factoring these rank-one matrices as row Gram matrices constructs an isometry from $\mathbb R^2$ to $\mathbb R^N$. Take its range for $P$ and put $F=\operatorname{diag}(1,\ldots,N)$. An orthogonal change puts either construction on the same prescribed mass block.

The compressed moments have equal scalar parts. Their traceless parts are the coefficients of $cB_b(z)/g(z)$ at $z^{-j-1}$. They first occur at $j=N-3$ with value $cZ$ in both examples. At $j=N-2$ they are $c(\sum_a a)Z$ and $c[(\sum_a a)Z+2X]$. Below order $2N-4$, every product of moments expanded in scalar and traceless terms has at most two traceless factors since $3(N-3)>2N-5$. The only possible trace pairings use two leading factors or one leading and one next factor, whose inner products agree. A word with at least one $P$ expands using $Q=I-P$ into such compressed products, with gaps summing to its order. Terms with only $I$ use the same full spectrum. At the cutoff the traces of $M_{N-2}^2$ differ by $8c^2$, so the portals give inequivalent response data. [Exact rational checks](ch39-low-rank-portal-words.py)
