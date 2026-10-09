import sympy as sp
t,c=sp.symbols('t c',positive=True)
E11=sp.Matrix([[1,0],[0,0]]);E12=sp.Matrix([[0,1],[0,0]])
x=t*E12; y=t*(c*E11+E12)
a=[(2*x+y)/3,(-x+y)/3,(-x-2*y)/3]; u=(x+y)/3; b=[u,-u]
I=sp.eye(2); Z=sp.Matrix(2,2,sp.symbols('z0:4'))
def P1(z):
    return (z-a[1])*(z-a[2])+(z-a[0])*(z-a[2])+(z-a[0])*(z-a[1])
diff=sp.simplify(P1(Z)-3*(Z-b[0])*(Z-b[1]))
print("factorization residual",diff)
st=lambda M:M.H
d=3; S=sum(a,sp.zeros(2))
print("sum a",S)
L1=sum((bk*st(bk) for bk in b),sp.zeros(2)); R1=S*st(S)/d**2+sp.Rational(d-2,d)*sum((aj*st(aj) for aj in a),sp.zeros(2))
L2=sum((st(bk)*bk for bk in b),sp.zeros(2)); R2=st(S)*S/d**2+sp.Rational(d-2,d)*sum((st(aj)*aj for aj in a),sp.zeros(2))
print("C2.1(i) R-L",sp.simplify(R1-L1)); print("C2.1(ii) R-L",sp.simplify(R2-L2))
Q1=sum(((bk*st(bk))**2 for bk in b),sp.zeros(2)); A1=sum((aj*st(aj) for aj in a),sp.zeros(2))
RQ1=sp.Rational(2,d**2)*A1**2+sp.Rational(d-4,d)*sum(((aj*st(aj))**2 for aj in a),sp.zeros(2))
Q2=sum(((st(bk)*bk)**2 for bk in b),sp.zeros(2)); A2=sum((st(aj)*aj for aj in a),sp.zeros(2))
RQ2=sp.Rational(2,d**2)*A2**2+sp.Rational(d-4,d)*sum(((st(aj)*aj)**2 for aj in a),sp.zeros(2))
print("C2.3(i) R-L",sp.factor(sp.simplify(RQ1-Q1))); print("C2.3(ii) R-L",sp.factor(sp.simplify(RQ2-Q2)))
# KT conjecture 2.4 (first inequality), with sum a = 0
s2=sum((aj*aj for aj in a),sp.zeros(2)) - S*S/d**2
KT1=sp.Rational(d-6,d)*sum(((aj*st(aj))**2 for aj in a),sp.zeros(2))+A1**2/d**2+s2*st(s2)/d**2 \
  +sp.Rational(2,d)*sum((aj*(aj+S/d)*st(aj+S/d)*st(aj) for aj in a),sp.zeros(2)) - sp.Rational(4,d**3)*sum((aj*S*st(S)*st(aj) for aj in a),sp.zeros(2))
KT2=sp.Rational(d-6,d)*sum(((st(aj)*aj)**2 for aj in a),sp.zeros(2))+A2**2/d**2+st(s2)*s2/d**2 \
  +sp.Rational(2,d)*sum((st(aj)*st(aj+S/d)*(aj+S/d)*aj for aj in a),sp.zeros(2)) - sp.Rational(4,d**3)*sum((st(aj)*st(S)*S*aj for aj in a),sp.zeros(2))
print("C2.4(i) R-L",sp.factor(sp.simplify(KT1-Q1)),[sp.factor(e) for e in (KT1-Q1).eigenvals()])
print("C2.4(ii) R-L",sp.factor(sp.simplify(KT2-Q2)),[sp.factor(e) for e in (KT2-Q2).eigenvals()])
for nm,M in [("2.1i",R1-L1),("2.1ii",R2-L2),("2.3i",RQ1-Q1),("2.3ii",RQ2-Q2)]:
    print(nm,"eigs",[sp.factor(e) for e in sp.simplify(M).eigenvals()])
