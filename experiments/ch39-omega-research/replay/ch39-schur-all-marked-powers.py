"""Separate the two small transfer branches even with the geometric endpoint metrics."""
from functools import lru_cache
from pathlib import Path
import runpy
from flint import arb,acb,acb_mat

work=Path(__file__).resolve().parent
rep=runpy.run_path(str(work/"ch39-schur-repeated-transfer.py"))
cov=rep["cov"]
model=rep["model"]
s=model["s"]
W=s["wall"]
J=s["inside"]
block=s["block"]
modes=cov["modes"]
square=rep["square"]
small=arb(".00365")
gap=arb(".04055")


@lru_cache(maxsize=512)
def sample(i,n):
    return acb(0,2*arb.pi()*i/n).exp()


def strict_positive(d,points=128):
    # Eliminate q. After squaring the remaining real polynomial has degree two.
    c=d[0,0].real
    a,b,e=d[1,0],d[0,1],d[1,-1]
    if not c>2*abs(a):
        return False
    h0=c*c+2*abs(a)**2-4*(abs(b)**2+abs(e)**2)
    h1=2*c*a-4*b*e
    h2=a*a
    allowance=(2*abs(h1)+8*abs(h2))*(2*arb.pi()/points)**2/8
    for i in range(points):
        z=sample(i,points)
        if not h0+2*(z*(h1+z*h2)).real>allowance:
            return False
    return True


def inner_symbol(m):
    pieces=[block(m,(7,6),J),acb_mat(2,len(J)),acb_mat(2,len(J))]
    for i,j,d,_ in s["shell"]["edges"]:
        if d in (724,804) and j in (7,6):
            row,col=(7,6).index(j),J.index(i)
            pieces[1+(d==804)][row,col]=pieces[0][row,col]
            pieces[0][row,col]=0
    inv=block(m,J,J).inv()
    poly={(0,0):-m[7,6]}
    for p,a in zip(modes,pieces):
        for q,b in zip(modes,pieces):
            key=p[0]-q[0],p[1]-q[1]
            value=(a*inv*b.conjugate().transpose())[0,1]
            poly[key]=poly.get(key,0)+value
    return {k:(v+poly[-k[0],-k[1]].conjugate())/2 for k,v in poly.items()}


def box(a,b,low):
    try:
        m=cov["matrix_u"](a.union(b))
        rows=model["fourier_rows"](m)
        base=square(rows[0])
        if not low:
            contrast=square([x-y for x,y in zip(rows[1],rows[2])])
            d={k:gap**2*v-contrast[k] for k,v in base.items()}
            return cov["positive"](d,128) or strict_positive(d)
        real_x=inner_symbol(m)
    except ZeroDivisionError:
        return False
    multiplier=arb(".059")-arb(".051")*(a+b)/2
    for row in rows[1:]:
        light=square(row)
        d={k:light[k]-small**2*v+multiplier*real_x[k] for k,v in base.items()}
        if not strict_positive(d):
            return False
    return True


def cover(a,b,low):
    if box(a,b,low):
        return 1
    assert b-a>arb(".0000001"),(a,b)
    mid=(a+b)/2
    return cover(a,mid,low)+cover(mid,b,low)


if __name__=="__main__":
    cut=arb(9)/10
    low=cover(arb(0),cut,True)
    high=cover(cut,arb(6400),False)
    q=arb(".165")
    distortion=(arb(7930000)/4730000)*arb(".04822").exp()/(1-q*q)
    assert distortion<arb("1.809")
    up_mid,up_low=arb(".00306")/arb(".8659"),arb(".00000629")/arb(".8573")
    up_light=(up_low*distortion).sqrt()
    up_split=(up_mid/distortion).sqrt()-up_light
    down=rep["required"](arb(".000265")/(arb(".01410")*distortion),
                         arb(".0000137")*distortion/arb(".01392"),2)
    assert up_light<small and up_split>gap and down>arb(".0565")
    print("spectral cells below and above 0.9",low,high)
    print("endpoint distortion, necessary up light ceiling and up split floor",
          distortion.str(12),up_light.str(12),up_split.str(12))
    print("two-pass down contrast floor",down.str(12))
