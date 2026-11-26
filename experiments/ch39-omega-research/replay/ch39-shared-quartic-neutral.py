"""A common scalar quartic fixes the charged spectrum and bounds a neutral pole."""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import arb, ctx, fmpq, fmpq_mat, fmpz_mat

root = Path(__file__).parent
graph = run_path(str(root/"ch39-shell-diagonal-gauge-links.py"))


def bisect(function,lo,hi):
    lo,hi = arb(lo),arb(hi)
    assert function(lo)<0<function(hi)
    for _ in range(82):
        mid=((lo+hi)/2).mid()
        value=function(mid)
        if value<0:
            lo=mid
        elif value>0:
            hi=mid
        elif hi-lo<arb("1e-20"):
            break
        else:
            raise ArithmeticError("Root enclosure is too wide to resolve the midpoint sign")
    return lo.union(hi)


def pole_polynomials(vertices,lengths,tree_count,*, edge_filter=None,gap_denominator=64):
    n=len(vertices)
    lap=fmpz_mat(n,n)
    for i,p in enumerate(vertices):
        for j in range(i):
            if graph["squared"](p,vertices[j]) in lengths and (edge_filter is None or edge_filter(p,vertices[j])):
                lap[i,j]=lap[j,i]=-1
                lap[i,i]+=1
                lap[j,j]+=1
    port=vertices.index((11,11,25))
    remain=[i for i in range(n) if i!=port]
    minor=fmpz_mat([[lap[i,j] for j in remain] for i in remain])
    assert minor.det()==tree_count
    for r in range(1,n+1):
        positive=fmpq_mat([[lap[i,j]+fmpq(1,n)-fmpq(int(i==j),gap_denominator)
                            for j in range(r)] for i in range(r)])
        assert positive.det()>0
    top,bottom=minor.charpoly(),lap.charpoly()
    common=top.gcd(bottom)
    return top//common,bottom//common


def pole_resolvent(vertices,lengths,tree_count,*, edge_filter=None):
    top,bottom=pole_polynomials(vertices,lengths,tree_count,edge_filter=edge_filter)
    def green(energy,unit,derivative=False):
        t=energy/unit
        if derivative == 2:
            p,r=top(t),bottom(t)
            dp,dr=top.derivative()(t),bottom.derivative()(t)
            return (p*r*bottom.derivative().derivative()(t)-top.derivative().derivative()(t)*r*r
                    +2*dp*r*dr-2*p*dr*dr)/(unit*r)**3
        if derivative:
            return (top(t)*bottom.derivative()(t)-top.derivative()(t)*bottom(t))/(unit*bottom(t))**2
        return -top(t)/(unit*bottom(t))
    return green


primary_logs = ("24.3145558253697","26.71304007023748","26.81725683017028")


def graph_constants(*, port_builder=pole_resolvent):
    layers=graph["site_choices"]()
    chosen=[min((a for a in layer if a[1:3]==target),key=lambda a:a[3])
            for layer,target in zip(layers,[(18,243),(18,27),(15,12)],strict=True)]
    x,_,common=graph["mass_solution"](chosen)
    unit=(2*x).exp()
    points=[p for p in product(range(30),repeat=3) if sum(y*y for y in p)==867]
    def orbit(p):
        return tuple(sorted(p))
    outer,middle,next_to,inner,pole=(1,5,29),(7,17,23),(11,11,25),(13,13,23),(17,17,17)
    colour=[p for p in points if orbit(p) in (outer,middle,next_to)]
    right=[p for p in points if orbit(p) in (middle,next_to,inner,pole)]+[pole]
    gc=port_builder(colour,(32,152,72,56),243)
    gr=port_builder(right,(56,12,68),12)
    return unit,common,gc,gr


def parameters(logs=primary_logs,*,include_graphs=False):
    ctx.prec=192
    data=run_path(str(root/"ch39-published-gauge-inputs.py"))
    log_q=data["q"](logs[2])
    shared=data["model"]["target"](11431,log_q)
    unit,common,gc,gr=graph_constants()
    g=(shared/15+common/(48*arb.pi()**2))**(-arb(1)/2)
    m3,m8=[(data["q"](s)-log_q).exp() for s in logs[:2]]
    d=(m8*m8-m3*m3).sqrt()
    b=((1-5*d.log())/21).exp()
    def a(h):
        return ((1+4*d.log()+h.log()+arb(15).log()/2)/20).exp()
    def mass_equation(h):
        return 15*a(h)**2-14*b*b-(d/h)**2
    h=bisect(mass_equation,".5","1.0")
    u=m3/(h*g)
    w=arb(28).sqrt()*b/g
    onsite=h*arb(15).sqrt()*a(h)
    vc,vr=1-21*a(h).log(),1-21*b.log()
    for v in (vc+4*d.log()+onsite.log(),(8*vc+6*vr+2*d.log()+8*onsite.log())/5):
        assert v.contains(0) and abs(v)<arb("1e-22")
    assert unit/64>m3*m3
    def active(energy,derivative=False):
        return g*g*(3*gc(energy,unit,derivative)/8+gr(energy,unit,derivative)/4)
    result = (g,u,w,m3,h,active)
    return (*result,unit,gc,gr) if include_graphs else result


def verify():
    g,u,w,m3,h,active=parameters()
    assert arb(".55888662")<h*h<arb(".55888663")
    def limiting(energy):
        return 1+2*w*w*(active(energy)+2/(9*u*u))
    bound=bisect(limiting,".0001",".003").sqrt()
    assert bound<m3
    print("lambda/g^2", (h*h).str(12),"singlet vacuum / Q",u.str(12),
          "neutral mass supremum / Q",bound.str(12),"/ triplet", (bound/m3).str(12))


if __name__=="__main__":
    verify()
