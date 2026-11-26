"""Full charged one-loop Pati-Salam thresholds with the Higgs fields at one shell site."""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

root = Path(__file__).parent


def vevs(log_octet):
    x = ((1+arb(15).log()/2+4*log_octet)/20).exp()
    w = ((1-5*log_octet)/21).exp()
    z = log_octet.exp()
    return x,w,z


def equation(t):
    x,w,z = vevs(t)
    return 15*x*x-14*w*w-z*z


def threshold(t):
    x,w,z = vevs(t)
    onsite = arb(15).sqrt()*x
    vcolour = 1-21*x.log()
    vright = 1-21*w.log()
    colour = vcolour+3*z.log()+z.log()+onsite.log()
    hypercharge = (8*vcolour+6*vright+2*z.log()+8*onsite.log())/5
    return hypercharge,colour


def verify():
    ctx.prec = 192
    # Positive contributions to the real Dynkin trace for a complex triplet
    # (including its conjugate) and the two charged right vectors.
    hyper = lambda dimension,y:2*dimension*Fraction(3,5)*y*y
    assert [hyper(3,Fraction(2,3)),hyper(3,Fraction(-1,3)),hyper(1,1)] == [
        Fraction(8,5),Fraction(2,5),Fraction(6,5)]
    graph = run_path(str(root/"ch39-shell-diagonal-gauge-links.py"))
    choices = graph["site_choices"]()
    trees = [(18,243),(18,27),(15,12)]
    graphs = [min((item for item in layer if item[1:3]==target),key=lambda item:item[3])
              for layer,target in zip(choices,trees,strict=True)]
    _,_,common = graph["mass_solution"](graphs)
    data = run_path(str(root/"ch39-published-gauge-inputs.py"))
    shared = data["model"]["target"](11431,data["q"]("26.81725683017028"))
    coupling = (shared/15+common/(48*arb.pi()**2))**(-arb(1)/2)
    for lo,hi,increasing in [("-.194","-.193",True),("1.715","1.716",False)]:
        low,high=arb(lo),arb(hi)
        sign=1 if increasing else -1
        assert sign*equation(low)<0<sign*equation(high)
        for _ in range(75):
            mid=((low+high)/2).mid()
            value=sign*equation(mid)
            if value>0:
                high=mid
            else:
                assert value<0
                low=mid
        t=low.union(high)
        x,w,z=vevs(t)
        assert equation(t).contains(0)
        assert all(term.contains(0) and abs(term)<arb("1e-22") for term in threshold(t))
        v_sigma=z/(4*coupling)
        v_right=arb(28).sqrt()*w/coupling
        assert ((4*coupling*v_sigma)**2+(coupling*v_right)**2/2-15*x*x).contains(0)
        print("effective vector determinant masses / Q",x.str(12),w.str(12),
              "physical octet / Q",z.str(12),"vevs / Q",v_sigma.str(12),v_right.str(12))


if __name__=="__main__":
    verify()
