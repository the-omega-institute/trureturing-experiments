"""Whole-tower relative determinant at zero energy on the 19-point shell."""
from functools import cache
from itertools import combinations, permutations, product
from flint import acb, acb_mat, arb, ctx
ctx.prec=192
R=867
outer=sorted(set(permutations((1,5,29))))
inner=sorted(set(permutations((13,13,23))))
sites=outer+inner+sorted(z for z in product(range(30),repeat=3) if sum(a*a for a in z)==R and z not in outer+inner)
@cache
def theta(d):
    return (1-arb(d)/(2*R)).acos()
def dot(a,b):
    return sum(x*y for x,y in zip(a,b))
def det(p,x,y):
    return sum(p[i]*(x[(i+1)%3]*y[(i+2)%3]-x[(i+2)%3]*y[(i+1)%3]) for i in range(3))
def unit(i,j):
    # Diagonal pole fixed under coordinate reflections.
    p=(17,17,17);x,y=sites[i],sites[j]
    z=acb(R*(R+dot(p,x)+dot(x,y)+dot(y,p)),-arb(R).sqrt()*det(p,x,y))
    return z/abs(z)
V=[20000*(sum(z)-35)*(sum(z)-49)**2-6-arb(2)*(sum(z)-35)/7 for z in sites]
edges=[]
for i,j in combinations(range(len(sites)),2):
    d=2*(R-dot(sites[i],sites[j]))
    if d<=804:
        edges.append((i,j,d,unit(i,j)*(-1 if d==200 else 1)))
assert len(edges)==135
rho=104*arb.pi()**2/(100*460**2*theta(804))
h0=acb.zeta_zero(35).imag
a,b=theta(724)-theta(244),theta(804)-theta(244)
def phases(h):
    return {724:acb(0,a*h).exp(),804:acb(0,b*h).exp()}
def matrix(u,ph):
    m=acb_mat(19,19)
    for j,w in enumerate(V):
        m[j,j]=u+w
    for i,j,d,g in edges:
        if u==0:
            diag=off=1
        else:
            v=(u*rho*theta(d)).sqrt()
            off=v/v.sinh()
            diag=off*v.cosh()
        m[i,i]+=diag;m[j,j]+=diag
        v=-off*g*ph.get(d,1)
        m[i,j]=v;m[j,i]=v.conjugate()
    return m
def logdet(u,ph):
    d=matrix(u,ph).det()
    assert d.imag.contains(0) and d.real>0
    return d.real.log()
if __name__=='__main__':
    reference=logdet(0,{})
    for n in (35,81150,116837,137352):
        h=acb.zeta_zero(n).imag
        diff=logdet(0,phases(h))-reference
        print('zero, relative log det, log-mass exposure shift',n,diff.str(16),(-diff/2).str(16))
    diff=logdet(0,phases(h0))-reference
    assert arb('.0166705899605351')<diff and diff<arb('.0166705899605352')
    for d in (arb(2)/15,arb(2)/3,arb(1)/3):
        print('relative one-loop inverse coupling', (d*diff/(4*arb.pi())).str(13))
