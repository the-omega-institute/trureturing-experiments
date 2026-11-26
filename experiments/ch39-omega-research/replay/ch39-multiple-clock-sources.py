"""Convex representation of all nonnegative clock stiffnesses and their scalar endpoint bounds."""
from itertools import combinations
from pathlib import Path
import runpy
from flint import acb,arb,acb_mat

work=Path(__file__).resolve().parent
rank=runpy.run_path(str(work/"ch39-clock-edge-rank-bound.py"))
geo=runpy.run_path(str(work/"ch39-symmetric-shell-actions.py"))
s,J,W=rank["shell"],rank["J"],rank["W"]
block=rank["block"]
m=s["matrix"](0,{})
inv=block(m,J,J).inv()
E=[J.index(i) for i in range(6)]
G=block(inv,E,E)
outer=(11,12,18)
A=block(m,outer,J)*block(inv,range(len(J)),E)
base=[rank["c"][0,0][W.index(i),W.index(j)] for i,j in zip(outer,outer[1:]+outer[:1])]


def base_radius():
    radii=[arb(0) for _ in range(3)]
    for count in range(1,7):
        for u in combinations(range(6),count):
            projected=block(A,range(3),u)
            update=projected*block(G,u,u).solve(projected.conjugate().transpose())
            for i in range(3):
                radii[i]=max(radii[i],abs(update[i,(i+1)%3]).upper())
    return radii


def norm_growth():
    maps=[block(m,W,J),acb_mat(7,len(J)),acb_mat(7,len(J))]
    for i,j,d,_ in s["edges"]:
        if d in (724,804):
            loc=W.index(j),J.index(i)
            maps[1+(d==804)][loc]=maps[0][loc]
            maps[0][loc]=0
    maps=[z*block(inv,range(len(J)),E) for z in maps]
    max_i=max((1+sum(sum(abs(z[W.index(j),i])**2 for j in (6,7,8)).sqrt() for z in maps)).upper()
              for i in range(6))
    max_o=max(sum(abs(A[k,i])**2 for k in range(3)).sqrt().upper() for i in range(6))
    return max_i*max_o


def robust_ratio(total,cond):
    tri=rank["tri"]
    # In addition to the original four inequalities, the outer-ring two subtriangles
    # must lie above minus the base bond angle.
    tests=tri["tests"]+[tri["terms"](tri["loop"](f),tri["y"]) for f in [(6,12,11),(7,12,11)]]
    uncertainty=total*norm_growth()
    todo=[(0,0,0)]
    count=0
    while todo:
        level,i,j=todo.pop()
        if any(tri["test_value"](t,level,i,j)<0 for t in tests):
            continue
        ph=[1]
        for k in (i,j):
            co,si=tri["trig"](level,2*k+1,1)
            ph.append(acb(co,-si))
        values=[abs(sum(x*y for x,y in zip(row,ph))) for row in rank["rows"]]
        lower=sorted(z.lower() for z in values)[1]-uncertainty
        upper=max(z.upper() for z in values)+uncertainty
        if lower > arb(".01947")*cond*upper:
            count+=1
            continue
        assert level<16,(level,i,j)
        todo.extend((level+1,2*i+a,2*j+b) for a in (0,1) for b in (0,1))
    return count


if __name__=="__main__":
    radius=base_radius()
    err=sum((r/abs(z)).asin() for r,z in zip(radius,base))
    phi=sum(v.arg() for v in base)+arb(0,err.abs_upper())
    assert 0<phi<arb(".10827")
    density=[]
    for f in geo["geometry"]["big"]:
        angle=phi if set(f)==set(outer) else (arb.pi()-phi).union(arb.pi())
        orig=geo["geometry"]["lifted"](f,geo["negative"])
        density.append((angle-orig)/geo["area"](f))
    diameter=arb("23.13440")
    span=max(v.upper() for v in density)-min(v.lower() for v in density)
    cond=arb(7930000)/4730000/(1-arb(".148")**2)*(diameter**2*span/2).exp()
    assert cond<arb("4.51")
    total=arb(".6")
    count=robust_ratio(total,cond)
    print("uniform base phase error, zero-mode condition product upper, cross norm per unit strength",
          err.upper().str(12),cond.upper().str(12),norm_growth().upper().str(12))
    print("admissible clock boxes for total stiffness up to",total.str(5),count)
