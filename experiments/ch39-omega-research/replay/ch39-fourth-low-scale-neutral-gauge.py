"""从 M_Z 到 10¹³ GeV 对逐能标自由荷六截止的两项中性规范残差作连续外包。"""
from pathlib import Path
from runpy import run_path

from flint import arb,ctx


def weak_bounds(root,p0,pmin,weights):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    axial,boson = lib("neutral-axial-top"),lib("neutral-bosonic")
    mp = arb("1.22089e19")
    rho = (arb("246.22")/mp)**2/2
    w = rho*arb(".63811")**2/2
    z = w+rho*3*arb(".467453")**2/10
    t = rho*arb(".8616")**2
    def entropy(mass):
        a = mass/pmin
        return a*(1+(1/(4*a)).log1p())
    pi = arb.pi()
    hi = ((2*axial["J"](w/pmin)+boson["j_neutral"](z/pmin))/(4*pi)
          +3*axial["J"](t/pmin)/(8*pi)+axial["L"](t/pmin)/(24*pi)+(entropy(w)+entropy(z))/(8*pi))
    coeff = (2*w+z)/(24*pi)+(w+z)/(48*pi)+t/(16*pi)+t/(24*pi)
    ref = arb("-.074256577194663365546").union(arb("-.074256577194663365544"))
    def high(x):
        return ref+arb(0,(hi+coeff/x).upper())
    low = lib("native-rh-ray-pvm")["weak_series"](56)
    def low_at(r):
        return low(r)+arb(0,(2*r*p0/pmin).upper())
    large = high(256*p0)
    whole = large
    for n in range(8):
        step = arb(2)**(n-6)
        for j in range(64):
            whole = whole.union(low_at(((64+j)*step).union((65+j)*step)))
    bound = sum(c*(large if i==15 else whole) for i,c in enumerate(weights)).upper()
    print("中能标的弱差高端余量、第16层为大截止时的弱差上界",hi.str(13),bound.str(13),flush=True)
    assert bound<0


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    weights,xmin,xmax,x18,p0,species,initial = lib("fourth-universal-neutral-gauge")["verify"](root)
    shared = lib("fourth-variable-scale-neutral-gauge")
    mp = arb("1.22089e19")
    pmin = (arb("1e6")/mp)**2
    weak_bounds(root,p0,pmin,weights)
    ctx.prec = 160
    edges = [xmin,4*p0,16*p0,64*p0,256*p0,arb("1e-27"),arb("1e-20"),arb("1e-15"),
             (arb("1.176397e14")/mp)**2,x18,xmax]
    gauge = lib("all-field-signed-crossing")
    rule = gauge["angular_rule"](7,20)
    factors = [(qq-arb(8)*col/3,mix-col) for _,_,(qq,mix,col) in species]
    def endpoint(p):
        return shared["angular_endpoint"](gauge,rule,species,p0,p,edges)
    previous = endpoint(p0)
    ranges = []
    for index in range(4,27):
        p = (arb(10)**(arb(index)/2)/mp)**2
        current = endpoint(p)
        low_end = len(edges)-1 if index<=12 else 4
        values = shared["residual_rectangle"](previous,current,weights,initial,factors,
                                              low_end=low_end,high_start=-1)
        ranges.append(values)
        print("上端能标的十进指数、两项残差连续外包",index/2,
              [(v.lower().str(10),v.upper().str(10)) for v in values],flush=True)
        previous = current
    return ranges


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
