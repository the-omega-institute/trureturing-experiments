"""Massive background-field moments including heavy scalars, vectors and Dirac fields."""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx, fmpq_mat

ROOT = Path(__file__).parent
MATCH = run_path(str(ROOT/"ch39-quartic-soft-adjoint-match.py"))
MODEL = MATCH["RATIO"]["MODEL"]
a = MATCH["a"]


def charged_blocks():
    tables = [MODEL["mass_blocks"](r) for r in (0,1,-1)]
    for label,(base,metric) in tables[0].items():
        weights = MODEL["weight_indices"](label[1])
        weight = sum(w*c for w,c in zip(weights,(5,3,-8)))
        if not weight:
            continue
        pos,neg = tables[1][label][0],tables[2][label][0]
        inv = fmpq_mat(metric.tolist()).inv()
        forms = [inv*fmpq_mat(m.tolist()) for m in (base,(pos-neg)//2,(pos+neg)//2-base)]
        rank = fmpq_mat(pos.tolist()).rank()
        yield a(weight),[arb_mat(m) for m in forms],pos.shape[0]-rank


def scalar_trace(r,blocks,sigma):
    q = r*r
    lam = a("1/12")
    ans = [arb(0)]*3
    for weight,(base,linear,quadratic),nullity in blocks:
        p = (base+r*linear+q*quadratic).charpoly()
        b,c,d = [p[nullity+j]/p[nullity] for j in (1,2,3)]
        power = -b,b*b-2*c,-b**3+3*b*c-3*d
        for j,x in enumerate(power,1):
            ans[j-1] += weight*x/lam**j
    for j in (1,2,3):
        ans[j-1] += (6*((lam*sigma*q)**(-j)-(lam*q)**(-j))
                     -24*((lam*(1+sigma*q))**(-j)-(lam*(1+q))**(-j))-6)
    return ans


def spin_weights(j):
    integral = lambda n: F(factorial(n)**2,factorial(2*n+1))
    scalar = (integral(j)-4*integral(j+1))/(16*integral(j+1))
    vector = 3*scalar-integral(j)/(2*integral(j+1))
    assert scalar==F(1,8*(j+1)) and vector==-2-F(5,8*(j+1))
    return a(scalar),a(vector)


def heavy_moments(r,z,w,gamma,blocks,sigma):
    trace,det = (3+r*r)*z+w,(2+r*r)*z*w
    discr = (trace*trace-4*det).sqrt()
    minus,plus = (trace-discr)/2,(trace+discr)/2
    heavy = 3*z+w
    d,h = [r*(a(k)*z).sqrt() for k in ("13/10","7/2")]
    ell = ((gamma/624)**2+a("1/4")).sqrt()
    D = d*acb(((ell+a("1/2"))/2).sqrt(),((ell-a("1/2"))/2).sqrt())
    C,E = D+d,D-d
    P,W = -h*D/d,-h*C/d
    scalars = scalar_trace(r,blocks,sigma)
    moments,ratios = [],[]
    for j in (1,2,3):
        fs,fv = spin_weights(j)
        vector = 18*(minus**(-j)+plus**(-j)-heavy**(-j))
        bosons = (fs*scalars[j-1]+fv*vector)/3
        fermion = sum(s*abs(m)**(-2*j) for s,m in [(1,W),(-1,P),(1,E),(-1,D)])
        moments.append(bosons+fermion)
        ratios.append(abs(fermion/bosons))
    return moments,ratios


def verify():
    ctx.prec = 180
    blocks = list(charged_blocks())
    results = []
    for n,kind,r,v,state in MATCH["matched_points"]():
        if kind!="center":
            continue
        z,w = 1/state[0],1/state[1]
        moments,ratios = heavy_moments(r,z,w,acb.zeta_zero(n).imag,blocks,a("1/256"))
        assert a("2.74e-5")<ratios[0]<a("2.76e-5")
        assert a("1.46e-9")<ratios[1]<a("1.48e-9")
        assert a("4.42e-14")<ratios[2]<a("4.45e-14")
        H = moments[1]/moments[0]**2
        results.append(H)
        print(n,"total second-to-first ratio",H.str(12),
              "fermion/background", [val.str(5) for val in ratios])
    assert all(u-v>a("2.23e-7") for u,v in zip(results,results[1:]))
    print("Consecutive total ratios",[ (v-u).str(8) for u,v in zip(results,results[1:])])


if __name__ == "__main__":
    verify()
