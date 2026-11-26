"""Integer normalizations for a single coherent inverse response on the whole gauge band."""
from math import factorial
from pathlib import Path
from runpy import run_path

from flint import arb, acb, ctx

finite = run_path(str(Path(__file__).with_name("ch39-finite-arithmetic-height.py")))
hardy, weight = finite["hardy"], finite["weight"]


def full_band():
    n = 4096
    order = 32
    mass = sum(1/arb(k).sqrt() for k in range(1,n))
    factor = (arb(11644)**2+(2*order)**2).sqrt()/n
    tail = arb(n)/10570 + arb(1)/2
    for k in range(1,order+1):
        tail += abs(arb.bernoulli(2*k))*factor**(2*k-1)/factorial(2*k)
    tail /= arb(n).sqrt()
    second = mass*arb(64).log()**2
    eps = arb("3.65e-22")
    witness = arb("11204.2")
    assert min(abs(hardy(witness)),abs(finite["finite_value"](witness))) > 17
    certified = 0
    for k in range(10570,11644):
        left, right = arb(k),arb(k+1)
        queue = [(left,right,abs(hardy(left)),abs(hardy(right)))]
        while queue:
            a,b,va,vb = queue.pop()
            # Linear interpolation of the phase-centered finite Dirichlet sum.
            if max(va,vb)+2*tail+2*eps+second*(b-a)**2/8 < 18:
                certified += 1
                continue
            assert b-a>arb(1)/1024
            mid = (a+b)/2
            v = abs(hardy(mid))
            queue.extend([(a,mid,va,v),(mid,b,v,vb)])
    print("Arithmetic magnitude exceeds 17 and stays below 18 on",certified,"height intervals")


def verify():
    ctx.prec,ctx.cap = 112,2
    full_band()
    root = acb.zeta_zero(10875).imag
    best = weight(root)
    eps, slack, near = arb("3.65e-22"),arb("1e-10"),arb(1)/1024
    rival = arb("10617.645332791202380").union(arb("10617.645332791202395"))
    t = abs(hardy(rival))+arb(0,eps)
    assert weight(rival)/best-1 > t/(13*21+12*t)+arb("5e-7")
    # Monotonicity in normalization reduces the entire feasible set to the worst endpoint.
    stack = [(arb(10603),(root-near).upper()),((root+near).lower(),arb("10621.28"))]
    while stack:
        a,b = stack.pop()
        g = a.union(b)
        gain = weight(g)/best
        if gain<1-slack:
            continue
        t = abs(hardy(g))+arb(0,eps)
        if gain < (1-slack)*(1+t/(13*20+12*t)):
            continue
        assert b-a>arb("1e-9")
        mid = (a+b)/2
        stack.extend([(a,mid),(mid,b)])
    c = 1/(13*20+12*arb(".02"))
    assert c>arb(".0038")
    ctx.prec = 240
    grid = 2**64
    center = int((root*grid).floor().unique_fmpz())
    scores = {alpha:{} for alpha in (18,19,20)}
    for m in range(center-12,center+14):
        g = arb(m)/grid
        t = abs(finite["finite_value"](g))
        for alpha, values in scores.items():
            values[m] = weight(g)/(1+t/(13*alpha+12*t))
    for values in scores.values():
        n = max(values,key=lambda m:values[m].mid())
        assert n==195923358475422178854846 and all(values[n]>v for m,v in values.items() if m!=n)
        l = values[n].lower()/best
        assert l>1-slack
        radius = (1-l+l*c*eps)/(5*l*c-arb(".015"))
        assert 0<radius<12/arb(grid)
    print("Admissible integer normalizations for one direct coherent call: 18, 19, 20")


if __name__ == "__main__":
    verify()
