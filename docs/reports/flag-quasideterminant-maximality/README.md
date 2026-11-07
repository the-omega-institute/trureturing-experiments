# Flag quasideterminant maximality

Exact checks for the maximality request in Section 9 item 4 of Dmitry V.
Talalaev, *Commuting Families of Flag Quasideterminants in Reflection Equation
Algebras of Type A*, arXiv:2610.32521v1. The coefficient field is Q(q), with q
indeterminate, and the algebra is the full reflection equation algebra, without
central-character relations. This report supports a preregistration outline;
it does not contain a Lean proof or claim a completed formal settlement.

Run `python3 check.py` with Python 3 and SymPy. The successful output is
`check.out` and ends in `ALL_CHECKS_OK`. Arithmetic is exact throughout.

The program constructs R_F, the flip P, Rcheck=R_F P and L tensor 1, and checks
that all sixteen reflection-equation entries and the six proposed relations
have the same row space over Q(q). It checks both directions of the
presentation obtained by eliminating d=q^(-2)(t-a), including all four
quadratic critical overlaps. The corresponding iterated Ore presentation has
normal order a,t,b,c and relations cb=bc+p, where
p=(1-q^(-2))((1+q^(-2))a^2-q^(-2)at). This is the standard polynomial-domain
Ore-extension route to the faithful fraction division ring; its construction
and the domain theorem remain obligations of a subsequent Lean implementation.

For fidelity, the checker reproduces source Theorem 2.1 in rank two with the
literal flags I1=b, I2=c-d*b^(-1)*a. Exact normal forms verify their commutativity,
b*I2=I2*b=Delta=bc-ad+(1-q^(-2))a^2, and centrality of t=a+q^2*d. The inverse
symbol B represents b^(-1); its rewrite rules follow from the Ore relations.

For independence, substitution a=0,b=u,c=v,d=w annihilates every original
relation and sends (b,Delta,t) to (u,uv,q^2*w). The exponent map
(i,j,k) -> (i+j,j,k) has the symbolic left inverse (r,s,t) -> (r-s,s,t).
For a nonzero F(X,Y,Z), choose N at least its Y-degree. Then
G=X^N*F(X,Y/X,Z) is a nonzero polynomial: its exponent map
(i,j,k) -> (N+i-j,j,k) has left inverse (r,s,t) -> (r+s-N,s,t).
Its specialization is u^N*F(u,v,q^2*w), again nonzero. Faithfulness therefore
shows that b,I2,t are algebraically independent in the fraction division ring.
Clearing the denominator of a hypothetical t in Q(q)(b,I2) contradicts that
independence. Thus the extra central element gives a proper commutative
extension in every faithful ambient localization containing the flags.

The symbolic inverse formulas apply to arbitrary exponent triples; bounded
monomial samples are additional checks. No extension of the substitution to
the noncommutative fraction division ring is used. The program also checks a
nonzero PBW commutator, including at q=2, to exclude a scalar or commutative
quotient as the proposed counterexample algebra.
