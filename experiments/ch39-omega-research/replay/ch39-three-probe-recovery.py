"""Global inverse bounds for a shared unknown kinetic scale and two independent spectral weights."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

root=Path(__file__).resolve().parent
pair=run_path(str(root/"ch39-arbitrary-probe-pair.py"))
s=pair["s"]
B,base,rounded,multiply=s["B"],s["base"],s["rounded"],s["multiply"]
LO,HI=Q(1,2048),Q(1,1152)
STEPS=(Q(1,128),Q(1024),Q(5120))


def kinetic_row(polys,poles,eta,c):
    x=B(eta)/c
    prefactor=multiply(B(pair["unit"])/c,*(x-B(j*j) for j in base["OFFSETS"] if j)) / multiply(*(x+v for v in poles))
    u=x/B(pair["unit"])
    out=[]
    for coeff in polys:
        v=B(0)
        for a in reversed(coeff):v=rounded(a+v*u)
        out.append(rounded(prefactor*v))
    return out


def norm(x):
    return max(-x.lo,x.hi)


def main():
    tri=s["ThreeFits"]()
    polys=[pair["numerator"](tri,k) for k in (0,1)]
    poles=[tri.poles[i] for i in tri.ids]
    small=Q(0)
    large=Q(0)
    # Each strip encloses the actual signed 2x2 response coefficients
    # throughout that scale interval, including at common endpoints.
    for j in range(320):
        c=B(LO+(HI-LO)*j/320,LO+(HI-LO)*(j+1)/320)
        low,second,third=[kinetic_row(polys,poles,eta,c) for eta in STEPS]
        small=max(small,sum(map(norm,low)))
        a,b=second
        d,e=third
        det=a*e-b*d
        assert det.lo>0
        large=max(large,max(norm(e)+norm(b),norm(d)+norm(a))/det.lo)
    assert small<Q("1e-58") and large<10**11,(float(small),float(large))
    mean=sum((w/tri.poles[i] for i,w in tri.fits[2].items()),B(0)).lo-Q("4e-55")
    assert mean/(HI+STEPS[0]/poles[0].lo)**2>1000
    assert 1/(4*LO*STEPS[1])<=Q(1,2)
    assert 1/(LO*poles[0].lo)<11
    residual_gain=48*2**219
    pre=base["data"]("ch39-cutoff-threshold-inverses.json")[1]
    assert Q(max(sum(map(abs,row)) for row in pre["inverse"]),2**pre["bits"])/(1-Q(1,2**100))<2**219
    l1=sum(norm(fit.get(i,B(0))-tri.fits[2].get(i,B(0))) for fit in tri.fits[:2] for i in tri.ids)
    assert l1<Q(".25")
    # S=scale difference, U=max two coefficients, E <= residual_gain*beta.
    # 1000*S <= eps + 11*E + small*U.
    # U <= large*(eps + E/1024 + S/2). Propagate the two input errors separately.
    for eps,E in ((Q(1),Q(0)),(Q(0),Q(1))):
        S=(eps+11*E+small*large*(eps+E/1024))/(1000-small*large/2)
        U=large*(eps+E/1024+S/2)
        assert S<Q(".00101")*eps+Q(".0111")*E
        assert l1*U+E<3*10**10*eps+2*10**8*E
    print("Three fixed probes recover c and the full 50-site weights throughout [1/2048,1/1152].")
    print("Upper bounds for the short row and long-pair inverse norm:",float(small),float(large))


if __name__=="__main__":
    main()
