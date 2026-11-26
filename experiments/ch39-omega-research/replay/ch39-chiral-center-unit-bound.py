"""任意逐模抬升中心偶质量时的共同质量单位上界。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

ROOT = Path(__file__).parent
prefix = run_path(str(ROOT / "ch39-monotone-chiral-phase-obstruction.py"))
source = prefix["source"]


def verify(capped_lift=False):
    model, _, _, initial, _ = prefix["torus"]["local"]["inputs"]()
    ctx.prec = 224
    bands = prefix["shell_bands"](model)
    initial[2] = 1 / arb(".117").union(arb(".119"))
    groups = [(q%2==0, "F", [2*x for x in w], bands[q][3:])
              for q,w in zip(model["charges"],model["weyl"])]
    groups += [(True, "V", [0,-7,-arb(21)/2], bands[0]),
               (False, "S", model["higgs"], bands[3][3:]),
               (False, "S", model["higgs"], bands[3][:3])]
    lmin,lmax = [(arb(mu)/model["mz"]).log() for mu in ("1e13","1e17")]
    max_lift = arb(8).sqrt()
    k0 = arb("2e13")
    ref = (k0/model["mz"]).log()
    directions = [(-115,333,-218),(-11,17,-6),(-7,-4,11),
                  (-5,3,2),(-8,5,3),(-7,12,-5),(-7,11,-4),(-4,11,-7),(-1,0,1),(0,1,-1),(-1,1,0),(0,-1,1),(1,-1,0),(1,0,-1)]
    kernels = {"S":model["scalar"],"F":lambda m,p:model["sheet"](m,p,False),
               "V":lambda m,p:model["sheet"](m,p,True)}
    dot = source["dot"]
    table = []
    for d in directions:
        assert sum(d)==0
        table.append((2*arb.pi()*dot(d,initial),dot(d,model["beta"]),
                      [dot(d,w) for _,_,w,_ in groups],arb.pi()*sum(map(abs,d))))
    def lower(a,b):
        # κ≥k0 迫使 L≥t+log(k0/M_Z)。当 κ 继续增加时，t 的左端无界。
        lower_L = lmin if capped_lift or a is None else lmin.max(ref+a)
        time = lower_L.union(lmax)
        p = (2*b).exp()
        low_p = (2*a).exp() if a is not None else None
        sharp, smooth = [], []
        for even,kind,_,levels in groups:
            sharp_sum = smooth_sum = arb(0)
            for h in levels:
                lo,hi = h.lower(),h.upper()
                upper_sharp = (b-lo.log()).max(0)
                upper_smooth = kernels[kind](lo**2,p)
                if (even and not capped_lift) or a is None:
                    lower_sharp = lower_smooth = arb(0)
                else:
                    top = max_lift*hi if even else hi
                    lower_sharp = (a-top.log()).max(0)
                    ref_ratio = ((2*(b-lmin)).exp()/lo**2).min(1) if capped_lift else (model["mz"]/(k0*lo))**2
                    lower_smooth = (kernels[kind](top**2,low_p)-ref_ratio/8).max(0)
                sharp_sum += lower_sharp.union(upper_sharp)
                smooth_sum += lower_smooth.union(upper_smooth)
            sharp.append(sharp_sum)
            smooth.append(smooth_sum)
        bounds = []
        for exp in (sharp,smooth):
            for n in (len(groups),len(groups)-1):
                bounds.append(max(((offset-slope*time-dot(w[:n],exp[:n]))/norm).lower()
                                  for offset,slope,w,norm in table))
        return bounds
    if capped_lift:
        # κh>M_Z 对每个中心奇模成立，含其最轻保留费米子。
        right = lmax + min(h.upper() for h in bands[1][3:]).log()
    else:
        right = (arb("1e17")/k0).log()
    pending = [(None,arb(-8)),(arb(-8),right.upper())]
    floor = [None]*4
    while pending:
        a,b = pending.pop()
        res = lower(a,b)
        if all(v>arb(".02" if capped_lift else ".15") for v in res):
            floor = [v if old is None else min(v,old) for old,v in zip(floor,res)]
            continue
        assert a is not None and b-a>arb("1e-6"), (a,b,res)
        mid = ((a+b)/2).mid()
        pending.extend(((a,mid),(mid,b)))
    print("偶模最多抬升 √8" if capped_lift else "κ≥2×10¹³ GeV、偶模任意抬升",
          "尖/有限减法下十九/十六标量的极差下界",
          *[v.str(9) for v in floor])


if __name__ == "__main__":
    verify()
    verify(capped_lift=True)
