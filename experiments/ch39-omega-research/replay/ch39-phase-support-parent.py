"""One-loop parent mass for arbitrarily small positive additions to the support."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

root = Path(__file__).resolve().parent
resolution = run_path(str(root/"ch39-phase-support-resolution.py"))
base, B = resolution["base"], resolution["B"]
unit, grid, secular = [base[k] for k in ("UNIT","grid","secular")]


def parent_endpoint(mass, ids, endpoints, kinetic_scale=None):
    p, modes, q2, target = base["physical_inputs"]()
    q = q2/B(mass**2)
    factor = B(1) if kinetic_scale is None else kinetic_scale
    poles = [grid(factor*c) for c in base["kinetic"](ids)]
    pure = sum(Q(hi,unit) < q.lo for _,hi in poles)
    assert 0 < pure < len(ids) and Q(poles[pure][0],unit) > q.hi
    assert q.hi < 8 and Q(poles[0][0],unit)*mass**2 > 1000**2
    assert all(h < l for (_,h),(l,_) in zip(poles,poles[1:]))
    residues = [[grid((B(1)-p)*w.get(i,B(0))) for i in ids] for w in endpoints]
    assert all(lo>0 or (lo==hi==0 and k>pure) for rho in residues for k,(lo,hi) in enumerate(rho))
    exposure = B(1)
    for lo,hi in poles[:pure]:
        exposure = (exposure*(q/B(Q(lo,unit),Q(hi,unit)))**25).positive_scale()
    active = 25*pure
    for mode, mult in zip(modes,base["MODES"]):
        for j in range(pure+1):
            brackets = []
            for rho in residues:
                args = (poles,rho,grid(p),grid(mode),grid(B(Q(125**2,mass**2))))
                lo = poles[j-1][1]+1 if j else 0
                hi = poles[j][0]-1
                for _ in range(90):
                    mid = (lo+hi)//2
                    a,b = secular(mid,*args)
                    if b < 0: lo = mid
                    elif a > 0: hi = mid
                    else: break
                brackets.append((lo,hi,args))
            lo = min(row[0] for row in brackets)
            hi = max(row[1] for row in brackets)
            # At fixed test energy the secular expression is affine in the weights.
            assert all(secular(lo,*args)[1] < 0 and secular(hi,*args)[0] > 0 for _,_,args in brackets)
            cell = B(Q(lo,unit),Q(hi,unit))
            assert cell.lo*mass**2 > 1000**2
            if cell.hi < q.lo:
                active += mult
                exposure = (exposure*(q/cell)**mult).positive_scale()
            else: assert cell.lo > q.hi
    ratio = exposure**11/target
    assert ratio.lo > 1 or ratio.hi < 1
    return (1 if ratio.lo > 1 else -1),active


def main():
    reference, moved, *_ = resolution["tables"]()
    for lo, weights, ids in [(105399287240169, [reference], sorted(reference)),
                             (124328988645806, [reference,moved], sorted(moved))]:
        assert parent_endpoint(lo,ids,weights)==(1,144)
        assert parent_endpoint(lo+1,ids,weights)==(-1,144)
        print(len(ids),"sites:",lo,"< parent mass / GeV <",lo+1,"with 144 open doublets")
    assert 124328988645806-105399287240170==18929701405636


if __name__=="__main__":
    main()
