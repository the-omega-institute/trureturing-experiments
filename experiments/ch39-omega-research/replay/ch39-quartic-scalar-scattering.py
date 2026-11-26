"""Scalar contact eigenvalue and matched thresholds for a common quartic coefficient."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import arb, ctx
from sympy import I, Rational, diag, eye, sqrt, zeros

ROOT = Path(__file__).parent
PLANCK = run_path(str(ROOT/"ch39-quartic-planck-branches.py"))
RATIO = PLANCK["RATIO"]
a = RATIO["a"]


def group_traces():
    # Use orthogonal Hermitian matrices and rational inverse norms.
    basis = [eye(5)]
    for j in range(1,5):
        basis.append(diag(*([1]*j+[-j]+[0]*(4-j))))
    for i,j in combinations(range(5),2):
        for phase in (1,I):
            h = zeros(5)
            h[i,j], h[j,i] = phase, phase.conjugate() if phase==I else 1
            basis.append(h)
    lift = RATIO["MODEL"]["LIFT"]
    norms = [(h*h).trace() for h in basis]
    casimir, wedge_casimir, adjoint = zeros(5), zeros(10), []
    for j,(h,norm) in enumerate(zip(basis,norms)):
        wedge = zeros(10)
        for pair,value in lift(h.todok()).items():
            wedge[pair] = value
        assert (wedge*wedge).trace() == (8 if j==0 else 3)*norm
        casimir += h*h/norm
        wedge_casimir += wedge*wedge/norm
        if j:
            adjoint.append((h,1/(2*norm)))
    assert casimir == 5*eye(5) and wedge_casimir == 8*eye(10)
    # Three real Gaussian pairings in E Tr Phi^4.
    moment = sum(w*z*(2*(h*h*k*k).trace()+(h*k*h*k).trace())
                 for h,w in adjoint for k,z in adjoint)
    assert moment == Rational(282,5)
    return moment


def scattering_certificate():
    moment = group_traces()
    sizes = [200,200,24,10]
    pure = [100*(100+1)+100*(10+10),100*(10+10),moment,5*(5+1)]
    assert [4*v/Rational(n) for v,n in zip(pure,sizes)] == [242,40,Rational(47,5),12]
    mixed = [4*1000, Rational(10*3*24,2), Rational(1,2)*10*4*5]
    assert [2*v/sqrt(sizes[0]*n) for v,n in zip(mixed,sizes[1:])] == [40,6*sqrt(3),2*sqrt(5)]
    def schur(p, lam):
        return p-242*lam-(40*lam)**2/(p-40*lam)-108*lam**2/(p-a("47/5")*lam)-20/(p-12*lam)
    cap = 8*arb.pi()
    assert schur(cap,a(".09728547906")) > 0 > schur(cap,a(".09728547907"))
    assert schur(a("21.76814"),a("1/12")) < 0 < schur(a("21.76815"),a("1/12"))
    print("Maximal tree contact coupling between 0.09728547906 and 0.09728547907")


def scalar_thresholds(products, r, lam):
    return [s-a(w)*lam.log()/12 for s,w in zip(RATIO["thresholds"](products,r),(90,92,93))]


def matches(lam):
    products = RATIO["product_polynomials"]()
    powers = PLANCK["mass_exponents"](products)
    match = RATIO["MODEL"]["matched"]
    for n,h,eps in RATIO["high_zero_readings"]():
        low,v,_ = match(scalar_thresholds(products,arb(1),lam),eps,arb(1))
        unit_low,unit_v,_ = match(RATIO["thresholds"](products,arb(1)),eps,arb(1))
        assert (v-unit_v).contains(0) and ((low/unit_low)**42*lam-1).contains(0)
        geom = (2*arb.pi()*(h*h+a("1/4")+30+a("2/15"))/(3*arb(5).sqrt()/2)).sqrt()
        def mass(r):
            return v/low*geom*PLANCK["mass_ratio"](r,powers)
        left,right = PLANCK["enclose_root"](mass,PLANCK["TARGET"],F(".4"),F(".5"),False)
        r = a(left).union(a(right))
        yield n,eps,r,scalar_thresholds(products,r,lam)


def verify():
    ctx.prec = 160
    scattering_certificate()
    for n,eps,r,scalar in matches(a("1/12")):
        low,v,lx = RATIO["MODEL"]["matched"](scalar,eps,arb(1),r)
        assert a(".4505") < r < a(".4506") and a(".774") < low < a(".776")
        assert low**2/arb.pi() < a("1/4") and a(".028") < lx.exp() < a(".032")
        assert v*eps*lx.exp() > 1000
        print(n,"u/v",r.str(9),"g_low",low.str(8),"v/GeV",v.str(8),"X",lx.exp().str(7))


if __name__ == "__main__":
    verify()
