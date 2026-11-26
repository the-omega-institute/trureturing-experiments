"""Kinetic and nonradial energy required by a common perturbative-strength quark fit."""
from pathlib import Path
from runpy import run_path
from flint import arb, acb, acb_mat

n = run_path(str(Path(__file__).with_name("ch39-shell-gw-native-arc-reference.py")))
m,arc = n["m"],n["arc"]
adj,hs,count,rows = n["adj"],n["hs"],n["count"],n["rows"]


def scalar_kinetic():
    s,k,_,_ = m["m"]["five"]["m"]["transported_spectra"]()
    coord = [(i,g) for i in range(6) for g in range(9)]
    def entry(i,g,j,b):
        if g==0 and b==0:
            return s[i,j]
        if g*b==0:
            return 0
        return k[8*i+g-1,8*j+b-1]
    p = m["c"]["p"]
    return acb_mat([[p[i].conjugate()*entry(i,g,j,b)*p[j] for j,b in coord] for i,g in coord])


if __name__=="__main__":
    t = n["relative_phase"](arb(26)**2)
    K = scalar_kinetic()
    output = arc["conjugate_output"](t)
    G,H = m["G"],adj(output)*output
    ratio = m["d2"]/m["u2"]
    x,b = arb(340),arb(23000)
    lower = arb("2.049")-x*m["u2"]/(4*arb.pi())**2
    assert count(rows(K+b*(H-ratio*G)-x*G),lower)==0
    target = arc["central_target"]()
    P = acb_mat(rows(m["F"])+rows(output))
    # Choose a weighted right inverse, preserving the same physical target.
    unit = m["m"]["eye"](54)
    Z = (K+arb("2.25")*unit).solve(adj(P))
    field = Z*(P*Z).solve(target)
    norm = hs(field,field).real
    rayleigh = hs(field,K*field).real/norm
    assert norm.sqrt() < 4*arb.pi() and rayleigh < arb("2.75")
    # Copy any one site's entire matrix to the graph in the native phase reference.
    phase = [acb(0,(-1)**i*t/2).exp() for i in range(6)]
    copy = acb_mat([[phase[i//9] if i%9==j else 0 for j in range(9)] for i in range(54)])
    assert count(rows(adj(copy)*K*copy/6),arb("1.347"))==9
    assert (arb("2.049")-arb("1.347"))/2 > arb(".3509")
    print("kinetic quotient above 2.049, constructed value",rayleigh.str(9),
          "with norm",norm.sqrt().str(9), "copied profile quotient below 1.347")
    assert (arb("2.049")-arb(".015908"))/2 > arb("1.0165")
