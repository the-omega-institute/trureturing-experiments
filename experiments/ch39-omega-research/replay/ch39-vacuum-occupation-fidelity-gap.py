"""第七壳真空命中所需的占据总变差、平方保真度损失与主本征线偏移。"""
from decimal import Decimal
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb


def check_distances(reference, alpha, derivative, sign, tv_radius, sine_radius, gap):
    def signed_interval(weight, left, right):
        lo, hi = derivative(weight,left,right)
        return (lo,hi) if sign > 0 else (-hi,-lo)

    upper, lower = [], []
    derivative_ranges = []
    center = arb(0)
    overlap = (1-sine_radius**2).sqrt()
    for p, weight in zip(reference,alpha,strict=True):
        low, high = signed_interval(weight,p.lower(),p.upper())
        center += p*(low+high)/2
        _, hi = signed_interval(weight,p.lower(),min(arb(1),(p+tv_radius).upper()))
        lo, _ = signed_interval(weight,max(arb(0),(p-tv_radius).lower()),p.upper())
        upper.append(hi)
        lower.append(lo)
        start = max(arb(0),(p.sqrt()*overlap-(1-p).sqrt()*sine_radius).lower())**2
        end = arb(1)
        if p < 1-sine_radius**2:
            end = (p.sqrt()*overlap+(1-p).sqrt()*sine_radius).upper()**2
        derivative_ranges.append(signed_interval(weight,min(start.lower(),p.lower()),
                                                  max(end.upper(),p.upper())))
    assert tv_radius*(max(upper)-min(lower)) < gap
    center = center.mid()
    variation = sum(p*max(abs(lo-center).upper(),abs(hi-center).upper())**2
                    for p,(lo,hi) in zip(reference,derivative_ranges,strict=True)).sqrt()
    spread = (max(hi for lo,hi in derivative_ranges)
              -sum(p*lo for p,(lo,hi) in zip(reference,derivative_ranges,strict=True)))
    assert 2*sine_radius*variation+sine_radius**2*spread < gap


def verify(root):
    def load(name):
        return run_path(str(root/("ch39-"+name+".py")))

    load("seventh-shell-vacuum-gradient-sign")["verify"]()
    witnesses = load("six-common-spectrum-moment-witnesses")["WITNESSES"]
    sieve = load("six-cyclic-line-sieve")
    adj, sq, row = (sieve[name] for name in ("adj","squared","row"))
    expint = load("chiral-minimal-neutral-gravity")["expint"]
    bounds = (
        ("0.03513294610",".043689",".037516"),
        ("0.03598788185",".023962",".035304"),
        ("0.03511179382",".016395",".017613"),
        ("0.01188942679",".006330",".008076"),
    )

    def check(root, levels, blocks, model, frames):
        a0 = levels[0]
        trace = sq(blocks[0])
        complement = sum(sq(m) for m in blocks[1:])
        gram = sum((adj(m)*m for m in blocks[1:]),acb_mat(2,2))
        leading = (1+(1-4*gram.det().real/complement**2).sqrt())/2
        assert arb(".99828") < leading < arb(".99830")
        fields = load("three-index-free-atoms")["broken_fields"](root,(2306888100092,))[2][0]
        epsilon = (arb(3000)/arb("1.22089e19")-a0).sqrt()
        shell, frame = model[0], frames[0]

        def vacuum(cutoff):
            x = cutoff**2
            return -sum(g*x*x*expint(3,z/x) for g,_,z in fields)/(4*arb.pi())

        def derivative(weight, left, right):
            lo = (a0+complement*left).lower()
            hi = (a0+complement*right).upper()
            assert 0 < lo <= hi
            exponential = sum(g*(-z/lo**2).exp().union((-z/hi**2).exp()) for g,_,z in fields)
            slope = -weight*complement*lo.union(hi)**3*exponential/(2*arb.pi())
            return max(arb(0),slope.lower()), slope.upper()

        for (labels,re,im), (gap_limit,tv,sine) in zip(witnesses,bounds,strict=True):
            rows = acb_mat([row(s["zero"],i-1).tolist()[0] for s,i in zip(model,labels,strict=True)])
            right = adj(rows)*(rows*adj(rows)).inv()
            null = acb_mat([[int(i==j) for j in range(7)] for i in range(7)])-right*rows
            pivot = max(range(7),key=lambda j:null[j,j].real.lower())
            direction = acb_mat([[null[j,pivot]] for j in range(7)])
            direction /= sq(direction).sqrt()
            cap = epsilon*sum(sq(acb_mat([[right[i,j]] for i in range(7)])).sqrt() for j in range(6))
            assert trace > cap**2
            source = acb_mat([[acb(x,y)] for x,y in zip(re,im,strict=True)])
            norm_squared = sum(x*x+y*y for x,y in zip(re,im,strict=True))
            value = arb(0)
            occupation, alpha = [], []
            for i, weight in enumerate(shell["w"]):
                if i in shell["triv"]:
                    r = row(shell["zero"],i)
                    h = abs((r*direction)[0,0])
                    error = epsilon*sum(abs(z) for z in (r*right).entries())
                    lo = a0+max(arb(0),(h*(trace-cap**2).sqrt()-error).lower())**2
                    hi = a0+(h*trace.sqrt()+error)**2
                    value += weight*vacuum(lo.lower()).union(vacuum(hi.upper()))
                else:
                    p = abs((acb_mat([frame[i]])*source)[0,0])**2/norm_squared
                    occupation.append(p)
                    alpha.append(-weight)
                    value += weight*vacuum(a0+complement*p)
            assert len(occupation)==12 and sum(occupation).contains(1)
            assert all(p>0 and a>0 for p,a in zip(occupation,alpha,strict=True))
            assert value>0 or value<0
            sign = 1 if value>0 else -1
            gap = (sign*(value-arb("1e-122"))).lower()
            assert gap > arb(gap_limit)
            tv_radius, sine_radius = arb(tv), arb(sine)
            check_distances(occupation,alpha,derivative,sign,tv_radius,sine_radius,gap)
            loss = max(Decimal(tv),Decimal(sine))**2
            assert loss >= Decimal("0.000065221776")
            print("低线、总变差下界、平方保真度损失下界", labels,tv,loss)
        # 对固定谱的主线，由两项分别估计两态的迹距离。
        assert ((arb(".006330")-1+leading)/leading)**2 > arb(".0000213")

    source = load("fifth-full-positive-joint")["verify"](root,return_operator="spectral")
    load("six-trivial-shared-endpoints")["experiment"](root,*source,check=check)


if __name__=="__main__":
    verify(Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent)
