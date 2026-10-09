# Krishna arXiv:2206.08154 Conjecture DUALSMALE (C*-algebraic Dubinin-Sugawa dual mean value) in A = C x C (sup norm).
# For P = prod (t - a_j) coordinatewise, critical points w satisfy P'(w) = 0 in both coordinates.
import sympy as sp, itertools
t=sp.symbols('t')
def check(a, z):
    n=len(a); polys=[sp.expand(sp.prod([t-ai[k] for ai in a])) for k in range(2)]
    der=[sp.diff(p,t) for p in polys]
    dz=[der[k].subs(t,z[k]) for k in range(2)]
    if all(d==0 for d in dz): return None
    normdz=max(abs(d) for d in dz); thr=sp.nsimplify(normdz/n)
    roots=[sp.roots(sp.Poly(der[k],t)) for k in range(2)]
    best=None; rows=[]
    for w0 in roots[0]:
        for w1 in roots[1]:
            w=(w0,w1); num=max(abs(polys[k].subs(t,z[k])-polys[k].subs(t,w[k])) for k in range(2)); den=max(abs(z[k]-w[k]) for k in range(2))
            q=sp.nsimplify(num/den); rows.append((w,q,bool(q>=thr)))
    return thr, rows
R=sp.Rational
for a,z in [([(0,R(3,2)),(3,R(3,2)),(3,R(3,2))],(0,0)),([(0,R(11,8)),(R(23,8),R(3,2)),(R(25,8),R(13,8))],(0,0))]:
    thr,rows=check(a,z); print("a=",a,"threshold ||P'(z)||/n =",thr)
    for w,q,ok in rows: print("   w=",w,"quotient=",q,float(q),"satisfies" if ok else "fails")
