# Exact check of a counterexample to Krishna's Higher Order C*-algebraic Smale Mean Value Conjecture (arXiv:2206.08154, label HIGHERMEAN)
# in A = C x C with the supremum norm, n = 3, z = 0, P = (x^3 + t x^2 + x, y^3 - 3y), t = 21 (and symbolic t > 18).
import sympy as sp
x,y,t=sp.symbols('x y t')
def check(tv):
    P1=x**3+tv*x**2+x; P2=y**3-3*y
    r1=sp.solve(P1,x); r2=sp.solve(P2,y)
    assert len(r1)==3 and len(r2)==3
    d1=[sp.diff(P1,x,k).subs(x,0) for k in range(4)]; d2=[sp.diff(P2,y,k).subs(y,0) for k in range(4)]
    nrm=lambda a,b: sp.Max(sp.Abs(a),sp.Abs(b))
    assert sp.simplify(nrm(d1[1],d2[1])-3)==0 and sp.simplify(nrm(d1[2],d2[2])-2*tv)==0
    cu=sp.solve(sp.diff(P1,x),x); cv=sp.solve(sp.diff(P2,y),y)
    worst=None
    for u in cu:
        for v in cv:
            val=nrm(P1.subs(x,0)-P1.subs(x,u), P2.subs(y,0)-P2.subs(y,v))
            lhs2=sp.nsimplify(nrm(d1[2],d2[2])/2)*val/9
            print(f'  w=({sp.nsimplify(u)}, {v})  ||P(z)-P(w)|| = {sp.N(val,12)}  k=2 lhs = {sp.N(lhs2,12)}  > 4: {bool(sp.N(lhs2)>4)}')
            assert sp.N(lhs2)>4
    return True
print('t = 21'); check(21)
print('symbolic: second coordinate gives |P2(0)-P2(v)| = 2 at v = +-1, so k=2 lhs >= (2t/2)*2/9 = 2t/9 > 4 iff t > 18')
print('OK')
