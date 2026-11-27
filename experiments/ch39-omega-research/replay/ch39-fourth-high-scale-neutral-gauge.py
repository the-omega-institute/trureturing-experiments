"""固定第四壳首枢轴，在高能半轴同时放开比较能标和十九个荷六截止。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    w,xmin,xmax,x18,p0,species,initial = lib("fourth-universal-neutral-gauge")["verify"](root)
    windows = lib("fourth-variable-scale-neutral-gauge")
    lib("neutral-residual-all-cutoffs")["verify"]()
    mp = arb("1.22089e19")
    pmin = (arb("1e17")/mp)**2
    # 与较低能标证书相同的高端对数矩误差仅用动量的下界，适用于任意更大动量。
    low_fixed = lib("native-rh-ray-pvm")["weak_series"](56)
    def low_uniform(ratio):
        return low_fixed(ratio)+arb(0,(2*ratio*p0/pmin).upper())
    high = lib("full-neutral-free-cutoff")["weak_column"]
    whole,large = lib("outer-shared-distance-kernel")["weak_bands"](low_uniform,high,p0,pmin)
    assert sum(c*(large if i==15 else whole) for i,c in enumerate(w)) < 0

    ctx.prec = 160
    edges = [xmin,4*p0,16*p0,64*p0,256*p0,arb("1e-27"),arb("1e-20"),arb("1e-15"),
             (arb("1.176397e14")/mp)**2,x18,xmax]
    assert all(a<b for a,b in zip(edges,edges[1:]))
    factors = [(qq-arb(8)*colour/3,mixed-colour) for _,_,(qq,mixed,colour) in species]
    gauge = lib("all-field-signed-crossing")
    rule = gauge["angular_rule"](7,20)

    def endpoint(p):
        return windows["angular_endpoint"](gauge,rule,species,p0,p,edges)

    def residual(left,right):
        return windows["residual_rectangle"](left,right,w,initial,factors,low_end=4,high_start=-1)

    previous = endpoint(pmin)
    certificate = []
    for k in range(35,43):
        # M = 10^(k/2) GeV，动量端点可用整数表示。
        p = arb(10)**k/mp**2
        current = endpoint(p)
        certificate.append(residual(previous,current))
        previous = current
    last = residual(previous,previous)
    # 对每个固定的容许截止向量界定从末端到任意更高的有限动量的变化。
    # 其第16、18层只需在待检动量满足约束，不要求沿途也满足弱差。
    variations = sum(abs(c) for c in w)
    tail = [12*xmax*variations*sum(abs(b[j]) for b in factors)/(4*arb.pi()*p) for j in range(2)]
    certificate.append([v+arb(0,t.upper()) for v,t in zip(last,tail,strict=True)])
    floors = [min(row[j].lower() for row in certificate) for j in range(2)]
    print("高能半轴两项规范残差下界、末端剩余变化上界",*[v.str(13) for v in floors+tail],flush=True)
    assert floors[0]>arb("769.0018") and floors[1]>arb("137.5984")
    assert tail[0]<arb(".550") and tail[1]<arb(".138")
    return floors


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
