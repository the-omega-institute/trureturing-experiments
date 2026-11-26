"""An interior positive solution for 22 consecutive arithmetic readings and all four shell phases."""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, arb_series, ctx

readings = run_path(str(Path(__file__).with_name("ch39-finite-nbb-coherent-four-windows.py")))["readings"]
indices=(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,22,34,38,66,67,75,82,352,357,363)
# These 26 actual heights give a positive basis.


def main():
    ctx.prec,ctx.cap=896,3
    gs=[v.imag for v in acb.zeta_zeros(1,400)]
    logs=[arb(2).log(),arb(3).log(),(arb(3)/2).log(),arb(6).log()]
    def col(g):
        return [g*g/(g*g+j*j) for j in range(22)]+[(2*g*t).cos() for t in logs]
    system=arb_mat([col(gs[i-1]) for i in indices]).transpose()
    weights=system.solve(arb_mat([[v] for v in readings(22)]+[[0]]*4)).entries()
    assert min(weights)>arb(".0014") and sum(weights).contains(1)
    # Use the actual unit-multiplicity spectrum on these first simple zeros.
    s=arb_series([arb(1)/2,1])
    xi=s*(s-1)/2*((-s/2)*arb.pi().log()).exp()*(s/2).gamma()*s.zeta()
    mass=xi[2]/xi[0]
    p=[1/(g*g*mass) for g in gs]
    tail=1-sum(p)
    lower=sum(v*(2*logs[-1]*g).cos() for v,g in zip(p,gs))-tail
    assert tail>0 and lower>arb(".18")
    print("Positive basis for the 22 consecutive readings and four doubled shell phases, minimum weight",min(weights).str(9))
    print("The complete canonical RH clock has its 2 log(6) cosine greater than .18; lower bound",lower.str(10))


if __name__=="__main__":
    main()
