"""Test every flavor-blind linear Higgs profile at the selected heights; fixed rank-one limiting profile."""
from itertools import combinations,product
from pathlib import Path
from runpy import run_path
from flint import arb,acb

m=run_path(str(Path(__file__).with_name("ch39-shell-octahedral-clock-modes.py")))


def coefficients(s,t):
    u=m["comps"](s,m["ends"][0])
    v=m["comps"](t,m["ends"][1])
    out=[]
    for j in range(3):
        a=[u[(j-k)%3]*m["zplus"]**k for k in range(3)]
        b=[v[(j-k)%3]*m["zminus"]**k for k in range(3)]
        c=[a[(k+1)%3]*b[(k+2)%3]-a[(k+2)%3]*b[(k+1)%3] for k in range(3)]
        norm=(abs(c[0])**2+abs(c[1])**2/2+abs(c[2])**2).sqrt()
        out.append([z/norm for z in c])
    return out


def hierarchy(n):
    h=acb.zeta_zero(n).imag
    s,t=[(h*f).sin()**2 for f in m["m"]["shells"][26][2]]
    cs=coefficients(s,t)
    target=arb(".061")
    total=0
    for fix in range(3):
        def cover(a,b,c,d):
            nonlocal total
            x=[a.union(b),c.union(d)]
            x.insert(fix,arb(1))
            A=[[abs(cs[j][(j-i)%3])*x[(j-i)%3] for j in range(3)] for i in range(3)]
            norm=sum(k**2 for row in A for k in row)
            lower=sum(abs(A[i][j]*A[k][l]-A[i][l]*A[k][j]).abs_lower()**2
                      for (i,k),(j,l) in product(combinations(range(3),2),repeat=2))
            if lower>3*target**2*norm**2:
                total+=1
                return
            assert b-a>arb(1)/2**14 or d-c>arb(1)/2**14
            if b-a>=d-c:
                mid=(a+b)/2
                cover(a,mid,c,d)
                cover(mid,b,c,d)
            else:
                mid=(c+d)/2
                cover(a,b,c,mid)
                cover(a,b,mid,d)
        cover(arb(0),arb(1),arb(0),arb(1))
    print("selected zero",n,"arbitrary linear profile second/top above .061; boxes",total)


def example():
    cs=coefficients(arb(0),arb(0))
    source=(1,arb(2).sqrt(),-1)
    M=[[cs[j][(j-i)%3]*source[(j-i)%3] for j in range(3)] for i in range(3)]
    assert abs(M[0][0])>0
    assert all((M[i][j]*M[k][l]-M[i][l]*M[k][j]).contains(0)
               for (i,k),(j,l) in product(combinations(range(3),2),repeat=2))
    # Exact vanishing follows from annihilating q_1 and q_2 in the common frame.
    print("profile moments (1,sqrt(2),-1); rank-one limiting entry",abs(M[0][0]).str(12))


if __name__=="__main__":
    example()
    for n in (25,255,719):
        hierarchy(n)
