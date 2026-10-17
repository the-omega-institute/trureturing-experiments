# Source algebra and signed interval inclusion

All identities here concern one exact coherent real source underlying
the saved interval arrays. Arbitrary independent endpoint selections are
not assumed to realize that source.

Let the symmetric arithmetic matrix Q have diagonal node matrix D and
two arithmetic generator columns G, with

\[
[D,Q]=\sigma GJG^T,\qquad J^T=-J.
\]

Partition the original indices into L, E and O. Assume the selected
principal E block A=Q_EE is invertible. Put

\[
\begin{aligned}
U&=A^{-1}G_E,&Y&=A^{-1}Q_{EL},\\
H&=Q_{LL}-Q_{LE}A^{-1}Q_{EL},&
R&=Q_{OL}-Q_{OE}Y,\\
g&=G_L-Q_{LE}U,&G'_O&=G_O-Q_{OE}U,\\
h&=R^TWG'_O.
\end{aligned}
\]

W is the original diagonal residual-weight matrix. Symmetry and
commutation with D_O suffice for the algebra. Positivity of W does not
turn the subtraction of its Gram into a positive contribution.

## Schur displacement identity

More generally let a symmetric current matrix M satisfy
`[D,M]=G Omega G^T` with skew Omega. Partition an invertible symmetric
principal block A from the remainder:

\[
M=\begin{pmatrix}A&B^T\\B&C\end{pmatrix},\qquad
L=BA^{-1},\qquad S=C-BA^{-1}B^T.
\]

The displacement block identities and inverse identity give

\[
\begin{aligned}
D_RB-BD_K&=G_R\Omega G_K^T,\\
D_KB^T-B^TD_R&=G_K\Omega G_R^T,\\
[D_K,A^{-1}]&=-A^{-1}G_K\Omega G_K^TA^{-1}.
\end{aligned}
\]

Expanding the commutator of the Schur subtraction therefore yields

\[
[D_R,S]=(G_R-LG_K)\Omega(G_R-LG_K)^T.
\]

This proof does not assume that D_K commutes with A; that omission would
be incorrect for a general 2x2 pivot. Applying the identity to the
arithmetic E block gives

\[
[D_L,H]=\sigma gJg^T,\qquad
D_OR-RD_L=\sigma G'_OJg^T.
\]

Taking the transpose of the second identity, and using the commutation
of W with D_O, gives

\[
[D_L,R^TWR]=\sigma(gJh^T+hJg^T).
\]

Scalar diagonal padding commutes with D_L. Thus the source

\[
C=(H-R^TWR)_{II}-(a_{\rm orig}+\rho a_{\rm reuse})I,
\qquad\rho=1/32,
\]

has four displacement generators `[g,h]_I` and skew kernel
`sigma[[J,-J],[-J,0]]`. The saved diagonal recovery adds back the complete
feature quadratic: the original far term and the joint boundary term
cancel algebraically. The original and reused padding scalars remain
distinct. An extra epsilon shift changes only the independent diagonal,
and multiplication by eta scales the skew kernel.

In detail, put V=Phi K Phi^T for the 112 original far features and
Q_Gamma=beta^T S_Gamma^-1 beta for the 19 boundary features. The fixed
source construction gives

\[
\begin{aligned}
\ell^{\rm orig}_n&=H_{nn}-(R^TWR)_{nn}-V_{nn}-a_{\rm orig},\\
\ell_n&=\ell^{\rm orig}_n+(Q_\Gamma)_{nn}
              -\rho(V_{nn}+a_{\rm reuse}),\\
FK_tF^T&=(1+\rho)V-Q_\Gamma.
\end{aligned}
\]

Adding the last quadratic diagonal to ell cancels both V and Q_Gamma,
leaving exactly the displayed diagonal of C. The displacement calculation
identifies every off-diagonal entry. This does not assume
`a_orig=a_reuse`, nor does it turn outward diagonal enclosures into
Loewner inequalities between arbitrary source realizations.

## Inertia and outward inclusion

For either sign of an invertible pivot block,

\[
\begin{pmatrix}I&0\\-BA^{-1}&I\end{pmatrix}
M
\begin{pmatrix}I&-A^{-1}B^T\\0&I\end{pmatrix}
=\operatorname{diag}(A,S).
\]

The congruence is invertible, so inertia adds across the pivot and
remaining Schur complement. A scalar pivot uses its strict sign. For a
real symmetric 2x2 block, negative determinant gives one sign each;
positive determinant and a strict trace sign give two equal signs.

An interval with encoded endpoints `[lo,hi]/2^b` that excludes zero has
reciprocal enclosure

\[
\left[\left\lfloor2^{2b}/hi\right\rfloor,
\left\lceil2^{2b}/lo\right\rceil\right]/2^b.
\]

This order is valid on either strict sign component. Integer floor and
ceiling, corner-product multiplication and exact outward grid lifting
preserve inclusion. At each accepted step the computed pivot inverse
contains the exact coherent inverse. Consequently the computed diagonal
and generator updates contain the exact Schur diagonal and generator
rows. The displacement identity restores the induction premise for all
remaining off-diagonal entries, without treating unrelated interval
endpoints as a coherent matrix.

Only complete removal of the full saved label set permits a final
inertia. Uncertain division or sign supplies no singularity conclusion
or final inertia from a prefix. For the saved eta=1 transcript every
scalar pivot interval is strictly positive, so the complete finite
congruence has inertia `(1263,0,0)`. The resulting source conclusion is
`C > epsilon I`. It does not provide an eta budget, the 131 feature and
two paid actions, the final bordered count or an infinite-dimensional
RH conclusion. These are ordinary algebraic explanations of the
interval computation, not compiled Lean proofs.
