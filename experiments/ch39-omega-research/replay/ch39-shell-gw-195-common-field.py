"""Test a single matrix Higgs and its covariant conjugate after the Wilson zero switch."""
from pathlib import Path
from runpy import run_path
from flint import acb, acb_mat, arb

folder=Path(__file__).resolve().parent
flavor=run_path(str(folder/"ch39-shell-gw-195-flavor.py"))
native=run_path(str(folder/"ch39-shell-gw-native-arc-reference.py"))
arc=native["arc"]
old=native["m"]
hs, adj, rows, count=(native[k] for k in ("hs","adj","rows","count"))


def operators():
    sing,family,vertex=flavor["scalar_models"]()
    phase=old["c"]["p"]
    coords=[(i,g) for i in range(6) for g in range(9)]
    def entry(i,g,j,b):
        if g==0 and b==0:
            return sing[i,j]
        if g*b==0:
            return 0
        return family[8*i+g-1,8*j+b-1]
    K=acb_mat([[phase[i].conjugate()*entry(i,g,j,b)*phase[j] for j,b in coords]
               for i,g in coords])
    F=acb_mat([[phase[i]*vertex[row,9*i+g] for i,g in coords] for row in range(9)])
    return K,F


def conjugate(F,t):
    phases=[acb(0,-(-1)**i*t).exp() for i in range(6)]
    return acb_mat([[F[row,j].conjugate()*phases[j//9] for j in range(54)] for row in range(9)])


if __name__=="__main__":
    t=native["relative_phase"](arb(676))
    K,F=operators()
    G=adj(F)*F
    target=arc["central_target"]()
    u2,d2=old["u2"],old["d2"]
    lower=(u2-62*d2)/arb("70.94")**2
    assert count(rows(G-62*G.conjugate()),lower)==54
    Y=conjugate(F,t)
    P=acb_mat(rows(F)+rows(Y))
    field=adj(P)*(P*adj(P)).solve(target)
    strength=hs(field,field).real.sqrt()
    assert strength<arb("11.593")<4*arb.pi()
    D=adj(Y)*Y
    bound=arb("2.0528")
    x,b=arb(350),arb(23000)
    assert count(rows(K+b*(D-d2/u2*G)-x*G),bound-x*u2/(4*arb.pi())**2)==0
    phases=[acb(0,(-1)**i*t/2).exp() for i in range(6)]
    copy=acb_mat([[phases[i//9] if i%9==j else 0 for j in range(9)] for i in range(54)])
    h=arb("1.3675")
    assert count(rows(adj(copy)*K*copy/6),h)==9
    assert (bound-h)/2>arb(".3426")
    unit=old["m"]["eye"](54)
    Z=(K+arb("2.4")*unit).solve(adj(P))
    trial=Z*(P*Z).solve(target)
    norm=hs(trial,trial).real
    assert norm.sqrt()<4*arb.pi() and hs(trial,K*trial).real/norm<arb("2.79")
    print("zero 195: ground-reference strength >70.94; arc-reference attainable",
          strength.str(12),"kinetic >2.0528; copied profile <1.3675")
