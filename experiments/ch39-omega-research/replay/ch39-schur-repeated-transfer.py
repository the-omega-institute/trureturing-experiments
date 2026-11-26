"""All powers of a fixed cyclic transfer constrained by its two lighter Fourier branches."""
from pathlib import Path
import runpy
from flint import arb

cov=runpy.run_path(str(Path(__file__).with_name("ch39-schur-ring-covariance.py")))
model=cov["shared"]
modes=cov["modes"]
c=arb(".0565")
weight=arb(2)/5


def square(row):
    result={}
    for p,v in zip(modes,row):
        for q,w in zip(modes,row):
            k=p[0]-q[0],p[1]-q[1]
            result[k]=result.get(k,0)+v*w.conjugate()
    return result


def box(a,b):
    try:
        rows=model["fourier_rows"](cov["matrix_u"](arb(a).union(b)))
    except ZeroDivisionError:
        return False
    base=square(rows[0])
    contrast=square([v-w for v,w in zip(rows[1],rows[2])])
    for row in rows[1:]:
        raw=square(row)
        expr={k:z-raw[k]/weight-contrast[k]/(c*c*(1-weight)) for k,z in base.items()}
        if not cov["positive"](expr,512):
            return False
    return True


def cover(a,b):
    if box(a,b):
        return 1
    assert b-a>arb(".000001")
    mid=(a+b)/2
    return cover(a,mid)+cover(mid,b)


def required(a,b,n):
    light=b**(arb(1)/n)
    return (a**(arb(1)/n)-light)/(1-light)


if __name__=="__main__":
    cells=cover(arb(0),arb(6400))
    m=model["shell"]["matrix"](arb(6400),{})
    eta=arb(12)/(6400+4-11)
    phi=(-m[6,12]).arg()
    rotation=2*(eta.asin()/2).sin()
    assert (1-eta)*abs(m[6,11])>abs(m[6,12])
    tail=2*arb(3).sqrt()*(phi.sin()+eta)/(3*phi.cos()-arb(3).sqrt()*phi.sin()-4*(eta+rotation))
    assert tail<c
    mu=arb(".165")/(1-arb(".165"))
    distortion=(1+2*mu)/(1-2*mu)
    ratios=[(arb(".00306")/arb(".8659"),arb(".00000629")/arb(".8573")),
            (arb(".000265")/arb(".01410"),arb(".0000137")/arb(".01392"))]
    for mid,small in ratios:
        bare=required(mid,small,2)
        framed=required(mid/distortion,small*distortion,3)
        print("required bare two-pass and framed three-pass contrasts",bare.str(12),framed.str(12))
        assert bare>c and framed>c
    print("spectral cells, compact contrast and tail upper",cells,c.str(5),tail.str(10))
