"""A compact-band Guinand-Weil calculation when a translated frequency overlaps zero."""
from pathlib import Path
from runpy import run_path
from flint import acb,arb,ctx

home=Path(__file__).resolve().parent
m=run_path(str(home/"ch39-shell-local-prime.py"))
geom=run_path(str(home/"ch39-shell-third-shell.py"))
ctx.prec=112
a,B=m["a"],m["B"]
pi=arb.pi()


def low(T,k):
    assert 2*a<k<B and k+B<arb(2).log()
    cs=(k*T).cos()
    C=B-k
    F0=m["density"](k,False)*cs
    # The two shifts have breakpoints where either sign crosses 0, 2a or 4a.
    knots=[arb(0)]+[x for off in (k,-k) for t in (-B,-2*a,arb(0),2*a,B) if (x:=t-off)>0]
    knots.sort(key=lambda x:float(x.mid()))
    assert all(x<y for x,y in zip(knots,knots[1:]))
    total=acb(0)
    for lo,hi in zip(knots,knots[1:]):
        mid=(lo+hi)/2
        active=[]
        for off in (k,-k):
            z=mid+off
            if abs(z)<B:
                active.append((off,1 if z>0 else -1,abs(z)<2*a))
        def integrand(u,analytic):
            f=acb(0)
            diff=acb(0)
            for off,sign,inner in active:
                d=m["density"](sign*(u+off),inner)
                f+=d*((u+off)*T).cos()/2
                if lo==0:
                    assert not inner
                    dk=(3*sign*C*C-3*C*u+sign*u*u)/(96*a**4)
                    diff+=(dk*cs+d*T*((off+u/2)*T).sin()*(u*T/2).sinc())/2
            if lo!=0:
                diff=(F0-f)/u
            return 4*f*(u/2).cosh()+diff*(u/2).exp()/(1j*u).sinc()
        total+=acb.integral(integrand,lo,hi,abs_tol=arb("1e-23"),rel_tol=arb("1e-23"))
    tail=(-knots[-1]/2).exp()
    total+=F0*((arb(1)/4).digamma()-pi.log()+2*(tail.atanh()+tail.atan()))
    assert total.imag.contains(0) and total.real.rad()<arb("1e-14")
    return total.real


def main():
    pairs=(((17,0),(29,1)),((17,1),(26,1)),((26,0),(29,0)))
    for n,edges in ((25,(".40384",".40465")),(255,(".38003",".38084")),(719,(".35640",".35824"))):
        T=acb.zeta_zero(n).imag
        D=m["denominator"](T)
        floor=arb(3)/8
        for (r,i),(s,j) in pairs:
            p,q=geom["shells"][r][2][i],geom["shells"][s][2][j]
            for f,v in ((4*p,1),(4*q,1),(2*abs(p-q),-2),(2*(p+q),-2)):
                if f<B:
                    value=low(T,f)
                else:
                    value,prime,hits=m["cosine"](T,f)
                    assert not hits and prime==0
                floor+=v*value/(16*D)
        assert arb(edges[0])<floor<arb(edges[1]) and floor.rad()<arb("1e-13")
        print("pole and archimedean terms, no primes, center",n,"normalized floor mean",floor.str(16))


if __name__=="__main__":
    main()
