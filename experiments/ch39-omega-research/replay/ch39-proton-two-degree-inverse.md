# Two contact probes at different degrees determine the overlap and scale

Fix a finite simple graph with common real hopping $t$. At every vertex use the same positive couplings $u,v$, with $w=v/u$ and $\rho=|c|^2$ both in $(0,1)$. Select sites of known degrees $z_1,z_2>0$ with $z_1\ne z_2$. Measure $A_j=X_{AB,j}$ at both, and $X_{AA,1},X_{BB,1}$ at the first. On any prescribed compact box inside the physical parameter domain these four raw measurements determine $(u,w,\rho)$ uniquely for all sufficiently small nonzero $t$. The threshold can depend on the box and sites. The nineteen-site shell has degrees $10,14,16,18$, giving six unordered degree pairs. In particular the extra mixed ratio distinguishes the two nearby branches of the [single-site fold](ch39-proton-shell-overlap-fold.md).

## Limiting measurements

Three combinations of the raw measurements suffice. Take $A_1$ together with
$
B_t=\frac{A_2-A_1}{(z_2-z_1)t^4},\qquad
C_t=\frac{X_{BB,1}-X_{AA,1}}{z_1t^8}.
$
Write $S(\rho)$ for the static mixed ratio. The fourth-order single-bond shift at $u=1$ is $\ell(w,\rho)S'(\rho)$, and $h(w,\rho)$ is the eighth-order coefficient of the single-bond difference $X_{BB}-X_{AA}$ at that scale. By the [coordination-number expansion](ch39-proton-overlap-jacobian-coordination-law.py), the three combinations converge together with their first parameter derivatives to
$
\Phi_0(u,w,\rho)=\big(S(\rho),\,S'(\rho)u^{-2}\ell(w,\rho),\,u^{-4}h(w,\rho)\big).
$
The static ratio and every shift along the static curve have equal $AA$ and $BB$ components. The eighth-order part of their difference is consequently $z_1$ times the single-bond value. All discarded terms after division tend uniformly to zero on the chosen box because the graph is finite and the static denominators are positive.

## Inverting the three combinations

Throughout $0<w,\rho<1$ the single-bond moments give $S'>0$, $h<0$ and
$
D=2h\,\partial_w\ell-\ell\,\partial_wh<0 .
$
For fixed $\rho$ the function $\ell/\sqrt{-h}$ thus has strictly positive derivative $-D/[2(-h)^{3/2}]$. Recover $\rho$ first from the static mixed ratio, then recover $w$ from $B_0/(S'\sqrt{-C_0})$. The scale is $u=(h/C_0)^{1/4}$ and $v=uw$. This inversion works even where $\ell$ vanishes. At $w=9/10$ it has opposite signs at $\rho=1/100$ and $\rho=1/2$.

For the sign checks, set $q=1+18\rho+45\rho^2+20\rho^3$ and $p=1+24\rho+90\rho^2+80\rho^3+15\rho^4$. After removing positive factors, the relevant polynomials are
$
H=\frac{49152w^4q^3p}{\rho(16-w^2)}\,h,\qquad
N=\frac{294912w^6q^3p^2}{\rho^2}\,D .
$
Their tensor Bernstein coefficient ranges on the unit square are $[-22127616,0]$ in bidegree $(2,10)$ and $[-17774007552,-768]$ in bidegree $(4,13)$. Every basis weight is positive in the open square. Both polynomials are therefore strictly negative there. The rational numerator and denominator of $S'$ have positive constant terms and nonnegative coefficients. The [exact certificate](ch39-proton-two-degree-signs.py) reconstructs all three functions from the six single-bond insertion moments. A [second check](ch39-proton-two-degree-inverse.py) uses the slice $\rho=1/2$ and counts the shell's degree pairs from its 135 edges.

## Small nonzero hopping and measurement precision

The explicit inverse has nonsingular derivative at every interior point. On a compact box $K$, injectivity and continuity of the derivative give a positive constant $m_K$ such that $|\Phi_0(x)-\Phi_0(y)|\ge m_K|x-y|$. For close pairs this follows by expanding about one endpoint, using the positive minimum singular value and uniform continuity of the derivative. For pairs at a fixed positive separation use compactness and injectivity. On the convex box, the Lipschitz constant of $\Phi_t-\Phi_0$ tends to zero by the derivative convergence. Hence $\Phi_t$ obeys the same bound with $m_K/2$ for all sufficiently small nonzero $t$. For the finite shell one can take the smallest threshold over the site pairs under consideration.

If each raw probability ratio has absolute error at most $\epsilon$, the errors in $A_1,B_t,C_t$ are bounded respectively by $\epsilon$, $2\epsilon/(|z_2-z_1||t|^4)$ and $2\epsilon/(z_1|t|^8)$. Together with the inverse bound this gives a parameter diameter of order $O_K(\epsilon |t|^{-8})$ for compatible points in the box. The degree difference and the common parameters across sites are used explicitly. This identifies the parameters of the finite contact probe under the specified measurement assumptions.
