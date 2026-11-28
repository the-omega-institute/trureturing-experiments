"""认证两项引力矩修复的局部解上，三支规范逆耦合的能标无关残差。"""
from math import factorial
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

ROOT = Path(__file__).resolve().parent
TARGETS = tuple(range(2_345_700, 2_345_747, 2))


def verify():
    repair = run_path(str(ROOT/"ch39-one-share-gravity-repair.py"))
    load = repair["load"]
    edges,kernel,moments,change,share = repair["setup"]()
    _, lib, spectra = load("two-zero-cutoff-gravity")["gravity_fields"](
        ROOT,(repair["SOURCE"],*TARGETS),tight_neutral=True)
    model,_,_,initial,_ = lib["chiral"]["source"]["mean"]["inputs"]()
    fields = [(0,0,0),*load("signed-gravity-unification")["gauge_weights"](model)]
    b1,b2,b3 = model["beta"]
    t = (b1-b3)/(b1-b2)
    direction = [1-t,t,-1]
    projection = arb_mat([[0,0,*direction]])
    baseline = sum(d*u for d,u in zip(direction,initial,strict=True))
    norm = abs(1-t)+abs(t)
    total = sum(abs(w) for _,_,w in edges)
    factor = 2/(4*arb.pi()*factorial(kernel["ORDER"]+1))
    gaps, predictions = [], []
    q = [sum(d*b for d,b in zip(direction,charges,strict=True)) for charges in fields]
    source = spectra[0]
    net = sum(w for _,_,w in edges)
    projected_source = projection*kernel["series"](source,fields,order=1)
    short_moments = arb_mat([[v] for v in [*moments.entries()[:4],moments.entries()[-1]]])
    first, relative, shifts = [], [], []
    for target,rows in zip(TARGETS,spectra[1:],strict=True):
        assert max(z.upper() for _,_,z in rows) < min(arb(1),arb(".0069")*min(a.lower() for a,_,_ in edges))
        series = kernel["series"](rows,fields)
        coeff = [arb_mat([[int(n==k) for k in range(5)]])*series for n in (0,1)]
        tails = [factor*sum(abs(2*a) for _,a,_ in rows),factor*sum(abs(g) for g,_,_ in rows)]
        args = (coeff,moments,change,[e*total for e in tails],[2*abs(edges[16][2])*e for e in tails])
        (s,d),_ = repair["repair_box"](args,share)
        gauge = repair["dilate"](projection*series,s)
        tail = factor*total*sum(sum(abs(v*b) for v,b in zip(direction,charges,strict=True)) for charges in fields)
        score = baseline+(gauge*(moments+d*change))[0,0]+arb(0,tail.upper())
        gap = abs(score)/norm
        # E1 的质量对数常数乘壳总权，同一分布对 x 的对数项在换重谱时抵消。
        prediction = net*sum(charge*(z/z0).log() for (_,_,z),(_,_,z0),charge in
                             zip(rows[1:],source[1:],q[1:],strict=True))/(4*arb.pi())
        # 一阶迁移式使用源见证矩，保留目标质量项与第一个逆矩。
        delta = projection*kernel["series"](rows,fields,order=1)-projected_source
        corrected = (delta*short_moments)[0,0]
        first.append(score*corrected>0)
        relative.append((abs(score-corrected)/abs(score)).upper())
        shifts.append((abs(score-prediction)/abs(score)).upper())
        print(target,"规范不变量、任意能标规范差下界、质量对数预报",score.str(9),
              gap.lower().str(9),prediction.str(9),corrected.str(9),flush=True)
        gaps.append(gap.lower())
        predictions.append(score*prediction>0)
    assert min(gaps)>arb(".01347") and sum(predictions)==8 and sum(first)==24
    assert max(relative)<arb(".2085")
    print("最小规范差下界、对数项与一阶式命中数、各自误差（一阶、对数）",min(gaps).str(9),sum(predictions),sum(first),max(relative).str(9),max(shifts).str(9),flush=True)


if __name__=="__main__":
    verify()
