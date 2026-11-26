"""Haar measure of the clock-edge orbit selected by the Gaussian determinant."""
from functools import lru_cache
from itertools import product
from pathlib import Path
import runpy
from flint import acb,arb,acb_mat

edge=runpy.run_path(str(Path(__file__).with_name("ch39-clock-edge-source.py")))
s=edge["s"]
pi=arb.pi()
representatives=[(0,7),(0,8),(1,6),(1,8)]
freq=[(a,b) for a,b in product(range(-3,4),repeat=2)
      if (a>0 or a==0 and b>0) and max(a,0)+max(b,0)<=3 and max(-a,0)+max(-b,0)<=3]
# Each plus-clock factor is in one of three inner columns, each minus in one of three inner rows.
N=7
roots=[acb(0,2*pi*i/N).exp() for i in range(N)]
normalization=s["matrix"](0,{}).det().real


def costs(i,j):
    m=s["matrix"](0,{724:roots[i],804:roots[j]})
    det=m.det().real/normalization
    g=m.inv()
    out=[]
    for v,w in representatives:
        phase=-m[v,w]
        r=(g[v,v]+g[w,w]-phase.conjugate()*g[v,w]-phase*g[w,v]).real
        out.append(det*r)
    return out


def coefficients():
    values={(i,j):costs(i,j) for i,j in product(range(N),repeat=2)}
    angles=list(values)
    tests=[]
    for l in (1,2,3):
        delta=[values[k][l]-values[k][0] for k in angles]
        const=sum(delta)/N**2
        terms=[]
        for a,b in freq:
            z=sum(v*roots[(-a*i-b*j)%N] for (i,j),v in zip(angles,delta))/N**2
            terms.append((a,b,2*z.real,-2*z.imag))
        tests.append((const,terms))
    return tests


@lru_cache(maxsize=131072)
def trig(l,c,r):
    n=1<<l
    angle=pi*arb(c%(2*n),r)/n
    return angle.cos(),angle.sin()


def enclosure(t,level,i,j):
    value,terms,hessian=t
    dp=dq=arb(0)
    for a,b,u,v in terms:
        co,si=trig(level,(2*i+1)*a+(2*j+1)*b,0)
        value+=u*co+v*si
        derivative=-u*si+v*co
        dp+=a*derivative
        dq+=b*derivative
    width=pi/(1<<level)
    error=width*(abs(dp)+abs(dq))+width**2*hessian/2
    return value+arb(0,error.abs_upper())


def density(tests,depth):
    tests=[(x,terms,sum((abs(a)+abs(b))**2*(u*u+v*v).abs_upper().sqrt() for a,b,u,v in terms))
           for x,terms in tests]
    yes=no=0
    cells=[(0,0,0)]
    while cells:
        level,i,j=cells.pop()
        inside=True
        for t in tests:
            val=enclosure(t,level,i,j)
            if val<0:
                no+=4**(depth-level)
                break
            if not val>0:
                inside=False
        else:
            if inside:
                yes+=4**(depth-level)
            elif level<depth:
                cells.extend((level+1,2*i+a,2*j+b) for a,b in product((0,1),repeat=2))
    return arb(yes)/4**depth,1-arb(no)/4**depth


if __name__=="__main__":
    tests=coefficients()
    # Check against a direct solve off the interpolation mesh, at the 68th height.
    p,q=edge["p"],edge["q"]
    raw,_,_,_=edge["evaluate"](p,q)
    m=s["matrix"](0,{724:acb(0,p).exp(),804:acb(0,q).exp()})
    for t,pair in zip(tests,representatives[1:]):
        val=t[0]+sum(u*(a*p+b*q).cos()+v*(a*p+b*q).sin() for a,b,u,v in t[1])
        assert (val-(raw[pair]-raw[0,7])*m.det().real/normalization).contains(0)
    lo,hi=density(tests,12)
    print("proportion selecting length-squared 724", (2*lo).str(12),(2*hi).str(12))
    assert 0<lo<hi<arb(".5")
