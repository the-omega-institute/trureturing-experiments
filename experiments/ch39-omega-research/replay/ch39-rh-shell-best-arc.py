"""Enclose the best observable 36-arc discrepancy under all 42 exact RH readings."""
from itertools import combinations, product
from pathlib import Path
from runpy import run_path

from flint import arb, acb, arb_mat, ctx

previous=run_path(str(Path(__file__).with_name("ch39-rh-shell-minimum-gap.py")))
directions=(((0,3),(1,2)),((0,2),(1,3)),((2,3),),((0,1),))
selected=(0,6,20,24)
signs=(-1,-1,1,1)


def dot(x,y):
    return sum(a*b for a,b in zip(x,y))


def raw(u,v):
    return [sum(u[i]*v[j]+u[j]*v[i] for i,j in places) for places in directions]


def shell_data():
    vec=[[1/arb(2).sqrt()]+[j/arb(34).sqrt() for j in p]
         for p in product(range(5),repeat=3) if sum(j*j for j in p)==17]
    return [(dot(u,v),raw(u,v),raw(u,u),raw(v,v)) for u,v in combinations(vec,2)]


def read_arc(row,z):
    target,a,b,c=row
    u,v=1+dot(b,z),1+dot(c,z)
    den=(u*v).sqrt()
    actual=(target+dot(a,z))/den
    grad=[x/den-actual*(y/u+w/v)/2 for x,y,w in zip(a,b,c)]
    return actual-target,grad


def main():
    ctx.prec,ctx.cap=1792,3
    gs=[acb.zeta_zero(i).imag for i in previous["support"]]
    lv=previous["prior"]["read"](42)
    logs=[2*arb(t).log() for t in (2,3,arb(3)/2,6)]
    design=arb_mat([[g*g/(j*j+g*g) for j in range(41)]+[(g*t).cos() for t in logs] for g in gs])
    dual=design.solve(arb_mat([[g*g/(41**2+g*g)] for g in gs])).entries()
    delta=dot(dual[:-4],lv[:41])-lv[-1]
    norm=sum(map(abs,dual[-4:]))
    n=[-b/norm for b in dual[-4:]]
    threshold=delta/norm
    rows=shell_data()
    arcs=[rows[i] for i in selected]
    # Find a finite trial direction, then place it exactly on the supporting
    # plane. The upper and lower bounds only use this final feasible point.
    z=[arb(x) for x in ("-.0000767","-.0000767",".0000268",".000249",".00001588")]
    for _ in range(4):
        values=[read_arc(r,z[:4]) for r in arcs]
        jac=arb_mat([df+[-s] for (_,df),s in zip(values,signs)]+[n+[0]])
        residual=[f-s*z[4] for (f,_),s in zip(values,signs)]+[dot(n,z[:4])-threshold]
        step=jac.solve(arb_mat([[v] for v in residual])).entries()
        z=[(a-b).mid() for a,b in zip(z,step)]
    rounded=[arb(int((t*2**140).floor().unique_fmpz()))/2**140 for t in z[:4]]
    point=[v*threshold/dot(n,rounded) for v in rounded]
    weights=design.transpose().solve(arb_mat([[f] for f in lv[:41]+point])).entries()
    assert min(weights)>arb(".0012") and sum(weights).contains(1)
    assert (sum(w*g*g/(41**2+g*g) for w,g in zip(weights,gs))-lv[-1]).contains(0)
    ceiling=max(abs(read_arc(r,point)[0]).upper() for r in rows)
    assert ceiling<arb(".000016")
    # Restrict all candidate competitors using the nine-point rigidity bound.
    radius=arb(1)
    for _ in range(5):
        radius=(148*ceiling*(1+3*radius)**2).upper()
    assert radius<arb(".0025")
    # For each competitor take the average Jacobian along the segment to
    # the trial point. Every entry lies in the interval matrix below.
    for _ in range(3):
        box=[arb(0,radius)]*4
        mat=arb_mat([read_arc(r,box)[1] for r in arcs])
        inv=mat.inv()
        # Start at zero to bound each coordinate by ||J^-1||_infty*ceiling.
        tighter=max(sum(abs(inv[i,j]) for j in range(4)).upper() for i in range(4))*ceiling
        radius=min(radius,tighter.upper())
    assert all(abs(z)<radius for z in point)
    mat=arb_mat([read_arc(r,[arb(0,radius)]*4)[1] for r in arcs])
    sensitivities=mat.transpose().solve(arb_mat([[v] for v in n])).entries()
    assert all(v*s>0 for v,s in zip(sensitivities,signs))
    floor=min((s*read_arc(r,point)[0]).lower() for s,r in zip(signs,arcs))
    assert floor>arb(".00001587918") and ceiling-floor<arb("1e-37")
    print("Optimum of the largest normalized shell overlap error",floor.union(ceiling).str(23))
    print("Constructive ceiling and global lower bound apart by", (ceiling-floor).str(5))
    print("At the feasible point the marked 14/17 shift is",read_arc(rows[4],point)[0].str(11))


if __name__=="__main__":
    main()
