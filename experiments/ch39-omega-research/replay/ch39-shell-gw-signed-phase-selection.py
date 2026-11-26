"""A positive and a negative log determinant on the fixed-angle CKM circle."""
from pathlib import Path
from runpy import run_path
from flint import arb

m=run_path(str(Path(__file__).with_name("ch39-shell-gw-mixed-determinant-cp.py")))
u,d=([x*x for x in m[key]] for key in ("up","down"))
t12,t13,t23=m["angles"]
beta=2*t12.sin()*t12.cos()*t23.sin()*t23.cos()*t13.sin()*(u[2]-u[1])*(d[1]-d[0])
s=u[0]+d[2]
# Coefficients of the central shifted characteristic polynomial.
A=m["determinant"](arb.pi()/2,0,1,1)
trace=sum(u)+sum(d)
B=m["determinant"](arb.pi()/2,1,1,1)-1-trace-A


def pole(z):
    return (((z+trace)*z+B)*z+A)/(beta*(z+s))


def root(c):
    lo,hi=arb(0),arb("1e-8")
    target=2*pole(0)-c
    assert pole(lo)<target<pole(hi)
    for _ in range(80):
        mid=(lo+hi)/2
        if pole(mid)<target:
            lo=mid
        else:
            assert pole(mid)>target
            hi=mid
    return lo.union(hi)


if __name__=="__main__":
    delta=arb("1.139")
    c=delta.cos()
    z=root(c)
    derivative=((3*z**2+2*trace*z+B)/beta-pole(z))/(z+s)
    sensitivity=z*derivative/delta.sin()
    curvature=delta.sin()**2/((pole(z)-c)*(pole(0)-c))
    assert arb("6.4664e-9")<z<arb("6.4665e-9") and sensitivity>70
    print("two CP branches for shifts between",root(arb(1)).str(10),"and",root(arb(-1)).str(10))
    print("central shift",z.str(11),"d phase/d log shift",sensitivity.str(11),"angular curvature",curvature.str(9))
