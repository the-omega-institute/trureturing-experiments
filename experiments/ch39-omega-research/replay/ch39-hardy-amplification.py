"""Sharp odd amplification counts at normalization 128 for both readouts."""
from math import factorial
from pathlib import Path
from runpy import run_path
from flint import acb, arb, ctx

model = run_path(str(Path(__file__).with_name("ch39-finite-arithmetic-height.py")))
hardy, weight = model["hardy"], model["weight"]


def penalty(t, degree, coherent):
    value = abs((degree*(t/128).asin()).sin())
    return value/(13+12*value) if coherent else value/208


def verify():
    ctx.prec, ctx.cap = 144, 2
    n, p = 4096, 32
    mass = sum(1/arb(k).sqrt() for k in range(1,n))
    bound = (arb(11644)**2+(2*p)**2).sqrt()
    tail = arb(n)/10570 + arb(1)/2
    for k in range(1,p+1):
        tail += abs(arb.bernoulli(2*k))*(bound/n)**(2*k-1)/factorial(2*k)
    component_norm = mass+tail/arb(n).sqrt()
    assert component_norm < 127
    assert mass > 126
    root = acb.zeta_zero(10875).imag
    best = weight(root)
    eps = arb("3.65e-22")
    near, slack = arb(1)/1024, arb("1e-10")
    witness = arb("10617.645332791202380").union(arb("10617.645332791202395"))
    z_witness = abs(hardy(witness))+arb(0,eps)
    excess = weight(witness)/best-1
    readouts = [(98,False),(99,False),(7,True)]
    for degree, coherent in readouts:
        assert degree*(z_witness/128).asin()<arb.pi()/2
        assert excess-penalty(z_witness,6 if coherent else 97,coherent)>arb(".000001")
        # The normalized response is bounded below in a full neighborhood of its zero.
        limit = arb(".02")
        lower = (degree/arb(128))*(1-(degree*(limit/128).asin())**2/6)
        effective = lower/(13+12*degree*(limit/128).asin()) if coherent else lower/208
        assert effective > arb(".00365")
    for j in range(-32,32):
        cell = root+near*(2*j+1)/64+arb(0,near/64)
        assert abs(hardy(cell))+eps < arb(".02")

    # A competitive point lies below the known height cap; the original tail bound persists.
    stack = [(arb(10603),(root-near).upper()),((root+near).lower(),arb("10621.28"))]
    while stack:
        lo,hi = stack.pop()
        cell = lo.union(hi)
        gain = weight(cell)/best
        if gain<1-slack:
            continue
        t = abs(hardy(cell))+arb(0,eps)
        if abs(t)<128 and all(gain < (1-slack)*(1+penalty(t,d,s)) for d,s in readouts):
            continue
        assert hi-lo > arb("1e-9")
        mid = (lo+hi)/2
        stack.extend([(lo,mid),(mid,hi)])

    ctx.prec = 240
    best = weight(root)
    scale = 2**64
    center = int((root*scale).floor().unique_fmpz())
    scores = {spec:{} for spec in readouts}
    for m in range(center-16,center+18):
        g = arb(m)/scale
        t = abs(model["finite_value"](g))
        for spec in readouts:
            scores[spec][m] = weight(g)/(1+penalty(t,*spec))
    for spec, rows in scores.items():
        n = max(rows,key=lambda m:rows[m].mid())
        assert n==195923358475422178854846 and all(rows[n]>v for m,v in rows.items() if m!=n)
        l = rows[n].lower()/best
        assert l > 1-slack
        c = arb(".00365")
        radius = (1-l+l*c*eps)/(5*l*c-arb(".015"))
        assert 0<radius<16/arb(scale)
    print("Signal normalization 128, finite component norm bound",component_norm.str(9))
    print("Standard odd degrees 99 and 7; canonical rotation queries 98 and 7; finite height numerator",center)


if __name__ == "__main__":
    verify()
