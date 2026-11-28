# Matching the complete fifth-order jet of three overlap cones in three steps

On the four-dimensional free torus $\prod_{j=1}^5 z_j=1$, $|z_j|=1$, let $\zeta=e^{2\pi i/5}$, $c=\cos(2\pi/5)$. The light centers have all $z_j=\zeta^s$ for $s=0,\pm1$ and the heavy centers have $s=\pm2$. Residual coordinates satisfy $\sum_j k_j=0$, with $z_j=\zeta^s e^{ik_j}$. One word step is one of the five translations or its inverse. A three-step real scalar and center-shift-invariant Clifford vector match all matrix coefficients through degree five, with a uniform free gap. In the subclass with invariant vector and a nondegenerate cone, three steps are minimal at degrees four and five.

## Scalar and vector before polar normalization

Put $M=\sum_j z_j/5$, $N=\sum_j z_j^2/5$, $t=2(1+c)$ and $f=-(t+6c)/(6t+28c)$. Take
$$
\begin{aligned}
 F={}&t-1+\Re\!\left(2tfM-t(1+f)M|M|^2-tfN\overline M-5fM^2+fN+(1+4f)M^3\right),\\
 B={}&-1+F+32(1-|M|^2),\qquad Y_j=\Im(z_j\overline M).
\end{aligned}
$$
Every scalar term takes at most three steps, and the components of $Y$ sum to zero. The vector is the very same function of $k$ around all centers, with linearization $k$. The values of $F$ at $s=0,1,2$ are respectively $0,0,5(1+\sqrt5)/4$ and inversion gives the remaining center values.

To check the Taylor coefficients in all directions, use permutation and inversion averages of the real Fourier characters. Frequencies $v\in\mathbb Z^5$ are equivalent modulo the all-ones vector. Subtracting the median coordinate minimizes their word norm $\sum|v_i|$. At a given center and along a direction $d$, the degree-$j$ moment for an orbit is the permutation average of $(v\cdot d)^j$ times $\cos(2\pi s\sum v_i/5)$ for even $j$ or $\sin(2\pi s\sum v_i/5)$ for odd $j$, up to a common nonzero factor at that degree. On the zero-sum hyperplane the symmetric polynomials in degrees two and three are spanned by $p_2,p_3$, where $p_j=\sum k_i^j$. Degree four is spanned by $p_4,p_2^2$. Evaluate on $(1,-1,0,0,0)$, $(1,1,-2,0,0)$ and $(1,1,-1,-1,0)$: $p_3$ is nonzero on the second and the first and third distinguish the quartic basis. Degree one vanishes by symmetry. The replay checks exact light-center equality in degrees $0,2,3,4$ over $\mathbb Q(\sqrt5)$, factoring out the common $\sin(2\pi/5)$ from odd rows. Inversion gives equality also at $s=-1$. Since $B<0$ at each light center, $Y/(-B)$ has the same coefficients through degree five: the first scalar discrepancy is $O(|k|^5)$ and $Y=O(|k|)$. The matrix obtained by polar normalization is analytic in this ratio there. [Exact orbit moments](ch39-three-step-fifth-jet.py)

## A gap on the entire torus

The sum of the absolute coefficients of the seven terms of $F$ is below seven. If $|M|\le4/5$ then $B>-1-7+32(1-16/25)>1$. If $|M|>4/5$, each phase lies within $\pi/2$ of $\arg M$ because its projection plus the other four unit projections must exceed four. If also $\|Y\|_2<1/200$, every such angular displacement $\delta_j$ satisfies $|\sin\delta_j|<1/160$, hence $|\delta_j|<1/100$. Subtract their mean to obtain $k_j$ with $\sum k_j=0$, $\max|k_j|<1/50$; the torus constraint makes their common remaining phase an aligned center. Every word of length at most three changes there by at most $3\max|k_j|$. Also $1-|M|^2\le\sum_j|e^{ik_j}-1|^2/5\le\max|k_j|^2$. These estimates bound the change of $B$ from its center value by $21/50+32/2500<1/2$. Its absolute center value is at least one. Altogether
$$
 \sqrt{B^2+\|Y\|_2^2}\ge1/200.
$$
To obtain linear normalization $k/m$ with $m>0$, use $Y/m$ and the bound $\min(1/2,1/(200m))$. A positive common rescaling of the vector and scalar preserves the normalized operator and can reach any prescribed bounded reference-gap threshold. A zero of the light operator requires $Y=0, B<0$. In the low-$|M|$ region $B>1$; in the high region $Y=0$ forces all five phases to coincide. The three indicated light centers thus exhaust all light zeros and have the same orientation. The two indicated heavy values have $Y=0,B>0$.

## The two-step obstruction and the boundary of the power construction

Suppose the Clifford vector is invariant under $z_j\mapsto\zeta z_j$ and has an invertible linearization at a light center. Near each light center the complete normalized matrix and $Y/(-B)$ determine each other analytically, since its normalized scalar component is nonzero. If the vector-ratio jets match to degree four, the scalar jets match to degree three: a first nonzero difference of the scalar reciprocals at degree $d\le3$ times the common linear vector produces a nonzero term of degree $d+1$. Permuting coordinates and inverting all phases preserves these linear scalar constraints, common negative light values and positive heavy values. Averaging therefore reduces the scalar test to real character orbits under permutations and inversion.

At radius two there are five such orbits. In the order-independent orbit basis the light-center difference rows of degrees zero, two and three span a space of rank three. Adding the difference of the heavy and light center values leaves that rank at three. This forces their averaged difference to vanish, contradicting the required signs. With radius three, the light differences through degree five have rank seven in the nine-dimensional orbit basis. The replay includes the two fifth-order directions $(1,1,1,1,-4)$ and $(2,1,0,0,-3)$. Adjoining the heavy contrast leaves the rank at seven. Thus a three-step scalar with equal fifth-order light jets cannot give opposite heavy signs. For invariant kinetics, sixth-order matrix matching consequently requires at least four steps. Four steps suffice by squaring the two-step selector vanishing through degree two at each light center and adding the invariant gap term as in the [finite-order construction](Ω复算说明.md#三个轻锥配平到任意指定有限阶需要多远的作用). The minimum in this subclass is therefore four at degree six.

The full space of radius-three averaged scalars has dimension nine. Its absolute moment rows at the light centers through degree four have rank nine. Consequently a three-step scalar vanishing to order five at all three light centers cannot retain strictly positive heavy values, even before averaging. This rules out obtaining the range bound $3\lceil K/5\rceil$ by raising such a flat scalar to powers. It leaves other high-order constructions open. For general two-step kernels the vector can change between centers, so their fourth-order feasibility also remains open. [Exact rank checks](ch39-three-step-fifth-jet.py)

## How exceptional would an unrestricted two-step solution be?

Even without averaging, a two-step scalar with equal third-order expansions at the three light centers is invariant under the center shift on the entire torus. Its 61 distinct frequencies have integer representatives $v$ with $\|v\|_1\le2$ (representatives differing by a nonzero multiple of the all-ones vector are at distance at least five). The 21 with $q=\sum_j v_j=0$ are invariant. Write $x=ik$ and $E_v=\exp(\sum_{j=1}^4(v_j-v_5)x_j)$, truncated to degree four. In the Taylor differences between $s=\pm1$ and zero, a frequency contributes $(\zeta^{sq}-1)E_v$. The rows of degrees zero through three have rank 40 modulo 31 at $\zeta=2$. This certifies a nonzero minor over the cyclotomic field, and the 21 invariant columns vanish identically. Equal scalar jets thus remove every charged mode.

For an arbitrary kinetic vector $C$, let $R$ be its central fourth-order Taylor polynomial in four independent output coordinates. A fourth-order match of normalized matrices requires $C_s-C_0=a_s R$ through degree four for $s=\pm1$, with $a_s$ the degree-three jet of $B_s/B_0-1$. Allowing these two amplitudes to vary freely gives a necessary linear system once $R$ is specified. Each of the four nonzero charges contributes coefficients summing to zero because $C$ vanishes at all five centers. For a charge $q$, choose any one frequency $v_q$ of that charge. Its independent difference columns are $(\zeta^{sq}-1)(E_v-E_{v_q})$ in each output component. This gives $4(40-4)=144$ columns. The amplitude columns are $-x^\beta R$ for $|\beta|\le3$ in each of two center blocks, another 70 columns. Reading the 69 positive-degree terms through order four in four outputs and two blocks gives 552 rows.

For $R$ from $Y_j=\Im(z_j\overline M)$ the rank is 214, modulo both 31 and 101. The replay uses $10iY$ to clear denominators, which rescales all amplitude columns by a single nonzero constant. Consequently there is a full-column-rank minor over characteristic zero, polynomial in the coefficients of $R$. Where it is nonzero, the only solution has zero charged vector coefficients. The vector is then invariant and the scalar sign obstruction applies. Any unrestricted solution must therefore lie in the proper algebraic set where every full-column-rank minor vanishes. This excludes a neighborhood of the displayed kinetic jet. Feasibility within the rank-deficient set remains open. [Two modular rank certificates](ch39-two-step-central-jets.py)

There is a second necessary condition directly on the scalar. Its values or its gradients must differ between some two light centers. To see this, write $b=-B$, and choose linear output coordinates so the linear term of the common ratio $C/b$ is $x$. Let $T_d$ map the 40 charged Fourier coefficients to the two Taylor differences through degree $d$. If $T_1 b=0$, cubic matching gives $T_3 C_j=[x_j T_3 b]_{\le3}$, since the scalar difference starts in degree two. These four identities and $T_1 b=0$ form a linear system in 200 unknowns. Its rank is at least 190 by modular elimination, and $T_2$ has rank 30 (the number of its rows). The full ten-dimensional kernel is therefore the subspace with all charged vector coefficients zero and $b\in\ker T_2$. For a fourth-order match the resulting invariant vector forces equality of the scalar third-order jets, and consequently equal signs at all five centers. Thus even arbitrary kinetic vectors cannot repair a two-step kernel with equal light-center scalar values and gradients. [Exact elimination](ch39-two-step-central-jets.py)

## Equal light-center Wilson masses force a fourth-order obstruction

Suppose that the three light values of $b=-B$ are equal. Fourth-order matching then fails for every two-step vector, including arbitrary anisotropic vectors and arbitrary scalar gradients. To prove this, identify the four-dimensional output with the tangent hyperplane $V=\{x\in\mathbb R^5:\sum x_j=0\}$ using the common invertible linear term of $C_s/b_s$. Write its common quotient jet as $x+F_2(x)+\cdots$. For any center-dependent scalar or vector field $G$ put
$$
 E(G)=(G_++G_-)/2-G_0,\qquad O(G)=(G_+-G_-)/2 .
$$
Let $a=E(b)_1$, $h=O(b)_1$ be linear forms. The degree-two vector contrasts are $E(C)_2=a x$ and $O(C)_2=h x$. At degree three, $E(C)_3\equiv aF_2$ and $O(C)_3\equiv hF_2$ modulo scalar polynomial multiples of $x$. Therefore $a O(C)_3-h E(C)_3\equiv0$.

Write $P$ for the zero-sum projection, $q=P(x_j^2)$, $u=P(x_j^3)$ and $L=\sum x_j^2/5$. For the two charged word families write
$$
 U_n=\sum_j u_j x_j^n,\qquad V_n=\sum_{i\le j}v_{ij}(x_i+x_j)^n .
$$
Here the complex vector coefficients represent the real part of the corresponding Fourier sum. Each family has coefficient sum zero because $C$ vanishes at all five uniform centers. Write $d_j=\cos(2\pi j/5)-1$ and $s_j=\sin(2\pi j/5)$. The Taylor contrasts are
$$
 E(C)_n=\Re\big[i^n(d_1U_n+d_2V_n)\big]/n!,\qquad
 O(C)_n=\Re\big[i^{n+1}(s_1U_n+s_2V_n)\big]/n!.
$$
In particular $d_1\Re U_2+d_2\Re V_2=-2a x$ and $s_1\Im U_2+s_2\Im V_2=2h x$. Polarize quadratic moments as symmetric bilinear forms. Expanding each word gives $U_3=U_1(u)$, $U_2(x,q)=U_1(u)-L U_1(x)$ and
$$
 V_3=3V_2(x,q)-2V_1(u)+6L V_1(x).
$$
Substituting the quadratic constraints into $E(C)_3,O(C)_3$ and using $(a x)(x,q)=(a(x)q+a(q)x)/2$ leaves
$$
 E(C)_3\equiv A u+L D x+\frac{d_2}{2s_2}h q,\qquad
 O(C)_3\equiv A' u+L D' x-\frac{s_2}{2d_2}a q
$$
for real linear maps $A,D,A',D'$ on $V$. Hence
$$
 H q\in W\quad\text{modulo multiples of }x,\qquad
 H=-\frac{s_2}{2d_2}a^2-\frac{d_2}{2s_2}h^2 ,
$$
where $W$ is spanned by $\ell(x) A u$ and $\ell(x)L D x$ with all linear forms $\ell$ and maps $A,D$.

The only quadratic forms $H$ satisfying this inclusion are the multiples of $L$. Wedge with $x$ and use the first four coordinates on $V$. The first family spans a space of dimension at most $64-4=60$ after wedging, since the four independent expressions $x u_j$ map to zero. Here independence refers to the four coefficient arrays giving $x u_j$. The second family is $L$ times arbitrary quadratic vectors. Its projected dimension is at most $40-4=36$ since the four linear multiples of $x$ vanish after wedging. Thus $\dim(x\wedge W)\le96$. Adding the ten columns for $x\wedge Hq$ increases the dimension by at most nine because $Lq$ already belongs to the second family. Evaluating on $\{-1,0,1,2\}^4$ gives ranks 96 and 105 modulo both 31 and 101. A nonzero minor lifts to characteristic zero and attains both upper bounds, so the kernel on quadratic forms is exactly the line spanned by $L$. [Exact projected ranks](ch39-two-step-central-jets.py)

Here $s_2>0$ and $d_2<0$. The two positive multiples of squares in $H$ give it rank at most two, whereas $L$ has rank four. Thus $H=0$ and both $a$ and $h$ vanish. The three scalar first jets now coincide, a case excluded by the rank-190 calculation above. Any fourth-order solution with unrestricted two-step kinetic vector must therefore have unequal Wilson masses at some pair of light centers. The three-step fifth-order construction has equal light masses, so three steps are minimal for fourth- and fifth-order matching under this equal-mass condition alone.

## At most twelve candidates on each line of fixed center masses

Fix an arbitrary real two-step scalar $b=-B$ with three positive light-center values $m_s=b_s(0)$ and two negative heavy-center values. No relation among the three positive values is required. For
$$
 g=\sin\theta_1+\cos\theta_1-\sin\theta_2-\cos\theta_2
$$
all five center values of $b+t g$ are independent of real $t$. On this entire line, at most twelve values of $t$ allow a two-step nondegenerate kinetic vector with fourth-order matching. The bound allows the kinetic vector to be chosen anew at each $t$.

A real two-step kinetic component $f$ has 61 coefficients. Let $V_m$ consist of those $f$ that vanish at the five centers and obey $m_0 Df_s=m_s Df_0$ at $s=\pm1$. The five value rows and twelve independent light-gradient rows have joint rank 17 in the unrestricted coefficient space by a finite-field minor. Hence the five value rows and eight gradient-ratio conditions have rank 13 for each choice of positive $m_s$, giving $\dim V_m=48$.

For any scalar $a$ define the linear map
$$
 K_a f=\left(J^4(a_0 f_s-a_s f_0)\right)_{s=-1,1},
$$
where $J^4$ truncates at degree four in the common local coordinates. On the 56-dimensional subspace of five-center-zero components, $\ker K_g$ is spanned by $g$. Indeed, the five value rows and the coefficient rows of $K_g$ together have rank 60 modulo 101 and 181, while $g$ lies in their kernel in characteristic zero. The Fourier replay specializes both a primitive fifth root and $i$ to each finite field. The gradient of $g$ at center $s$ is $(\cos(2\pi s/5)-\sin(2\pi s/5))(e_1-e_2)$. This factor is positive at zero and negative at $+1$. Since $m_0,m_1>0$, $g\notin V_m$, so $K_g$ restricted to $V_m$ is injective. [Exact directional rank and center interpolation](ch39-two-step-central-jets.py)

Choose 48 rows detecting this injection in a basis of $V_m$. The determinant of the same rows for $K_b+t K_g$ is a degree-48 polynomial with nonzero leading coefficient. Fourth-order matching with invertible velocity requires four independent kinetic components in its kernel. At each feasible $t_0$ the square submatrix thus has nullity at least four. A constant change of column basis makes four columns vanish at $t_0$ and proves that the determinant is divisible by $(t-t_0)^4$. There can be at most twelve distinct such parameters. In every affine space with the five masses prescribed, the candidate set is therefore contained in a proper algebraic subset of measure zero. Determinant roots still require four independent slopes and the gap and sign conditions on the full torus.
