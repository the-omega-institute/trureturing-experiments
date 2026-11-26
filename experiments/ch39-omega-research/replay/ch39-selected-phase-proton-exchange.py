"""26 路七三角相位在电荷五树图交换中允许抵消的最小反向份额。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx

root = Path(__file__).parent
heavy = run_path(str(root/"ch39-shell-current-word-floor.py"))
selector = run_path(str(root/"ch39-port-coupling-competition.py"))
G = heavy["green"]


def attaining_weights(phase):
    p,q = arb(5)/26,arb(21)/26
    rays = [G[1,13],G[1,17],phase*G[1,16],phase.conjugate()*G[1,1]]
    system = arb_mat([[1,1,0,0],[0,0,1,1],
                      [v.real for v in rays],[v.imag for v in rays]])
    weights = system.solve(arb_mat([[p],[q],[0],[0]])).entries()
    assert all(t>0 for t in weights)
    return weights


def verify():
    ctx.prec = 192
    selected = selector["prepare"]()
    phase = (selected[1,0]/abs(selected[1,0]))**5
    p,q = arb(5)/26,arb(21)/26
    n = G.nrows()
    # The cone also works for arbitrary arm masses.
    smallest = min(G[i,j].real.lower() for i in range(n) for j in range(n))
    spread = 2*arb(".12776994016814").atan()
    span = 2*arb.pi()-5*arb(".892")+spread
    assert smallest*(span/2).cos() > arb(".01120")

    weights = attaining_weights(phase)
    threshold = weights[3]
    assert arb(".13675218914066") < threshold < arb(".13675218914067")
    critical = {1,2,5}
    rest = set(range(n))-critical
    rotation = [heavy["points"].index((y,z,x)) for x,y,z in heavy["points"]]
    assert {1,rotation[1],rotation[rotation[1]]} == critical

    def project(form,z):
        return (form.conjugate()*z).real

    def support(row,form,r):
        minima = [min(project(form,u*z).lower() for z in row)
                  for u in (1,phase,phase.conjugate())]
        return p*minima[0]+(q-r)*minima[1]+r*minima[2]

    rows = [[G[i,j] for j in range(n)] for i in range(n)]
    special = acb(0,1)*(G[1,13]-G[1,17])
    special /= abs(special)
    generic = acb(-1,-29)
    generic /= abs(generic)
    for group,opposite,form in ((critical,rest,special),(rest,critical,generic)):
        for i in group:
            assert support(rows[i],form,arb(0)) > arb(".016")
            if group == rest:
                assert support(rows[i],form,threshold) > arb(".00013")
            for j in opposite:
                midpoint = [(a+b)/2 for a,b in zip(rows[i],rows[j])]
                assert support(midpoint,form,arb(0)) > arb(".019")
                assert support(midpoint,form,threshold) > arb(".0049")
    # On each critical pure-source row the four active rays saturate the
    # separating functional at the exact solution. Neutral equality is
    # perpendicularity to G[1,13]-G[1,17]; rotation supplies the other rows.
    for u,active in ((1,[13,17]),(phase,[16]),(phase.conjugate(),[1])):
        value = project(special,u*G[1,active[0]])
        assert all(project(special,u*G[1,j])>value
                   for j in range(n) if j not in active)

    print("Exact minimum reverse current weight",threshold.str(18),
          "fraction of the phased arm",(threshold/q).str(18))
    print("Attaining neutral, neutral, forward, reverse weights",*[v.str(15) for v in weights])


if __name__=="__main__":
    verify()
