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

The first minimum is strictly the smallest. In that cyclic order its edge indices are $(4,3,2,1)$, with cycle product $-\sigma$. Its reversed orientation adds the same real weight. All other contributions to the difference decay strictly faster, so the connected four-source responses differ for all sufficiently large $L$. This uses access to four separately located insertions in the canonical continuum probe. The static spectrum, bubble and scalar triangle already specified sample lower-order information. Distinguishing this pair by a fourth-order response leaves the identification of general portals a separate question.

[Exact matrix and loop checks](ch39-four-mass-portal-cut.py) · [Large-argument expansion of the modified Bessel function, DLMF 10.40.2](https://dlmf.nist.gov/10.40.E2)
