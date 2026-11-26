"""第 11661 枚三代边界在中心偶谱脱耦时的逐阈值两圈交点。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx

ROOT = Path(__file__).parent
chiral = run_path(str(ROOT / "ch39-chiral-common-gap-crossing.py"))
sm = run_path(str(ROOT / "ch39-two-loop-shell-separation.py"))
flow = run_path(str(ROOT / "ch39-two-tier-dirac-running.py"))["integrate_segments"]


def thresholds(model, layers, full):
    knots = []
    for q in (1,3):
        color = 3 if q==1 else 1
        indices = [arb(2*color*q*q)/60,arb(color)/2,arb(int(color==3))]
        casimirs = [arb(q*q)/60,arb(3)/4,arb(4)/3 if color==3 else 0]
        for j,h in enumerate(layers[q]):
            fermion = int(j>=3)
            scalar = int(q==3 and (full or j>=3))
            if not fermion and not scalar:
                continue
            one = [v*(4*fermion+scalar)/3 for v in indices]
            rep = model["weyl"][model["charges"].index(q)]
            expected = [2*fermion*w+scalar*a for w,a in zip(rep,model["higgs"])]
            assert all((a-b).contains(0) for a,b in zip(one,expected))
            two = [[indices[i]*(4*(fermion+scalar)*casimirs[k]+
                     ((20*fermion+2*scalar)*arb((0,2,3)[i])/3 if i==k else 0))
                    for k in range(3)] for i in range(3)]
            knots.append((h.log(),one,two,fermion,scalar))
    knots.sort(key=lambda row:float(row[0]))
    assert all(a[0]<b[0] for a,b in zip(knots,knots[1:]))
    return knots


def crossing(model,layers,initial,full):
    knots = thresholds(model,layers,full)
    if full:
        L = arb("28.70").union(arb("28.85"))
        t = arb("4.265").union(arb("4.29"))
    else:
        L = arb("28.99").union(arb("29.05"))
        t = arb("4.379").union(arb("4.385"))
    active = [row for row in knots if row[0]<t]
    assert all(row[0]<t or row[0]>t for row in knots)
    beta, B = sm["BETA"], sm["B"]
    phases = [(beta,B)]
    for _,w,add,_,_ in active:
        beta = [a+b for a,b in zip(beta,w)]
        B = [[x+y for x,y in zip(r,s)] for r,s in zip(B,add)]
        phases.append((beta,B))
    total = [sum(row[1][i] for row in active) for i in range(3)]
    pi = arb.pi()
    matrix = arb_mat([[sm["BETA"][i],total[i],2*pi] for i in range(3)])
    base = [2*pi*a+sum(row[0]*row[1][i] for row in active) for i,a in enumerate(initial)]
    def image(L,t,steps):
        durations = [L-t+active[0][0]]
        durations.extend(b[0]-a[0] for a,b in zip(active,active[1:]))
        durations.append(t-active[-1][0])
        assert all(v>0 for v in durations)
        segments = [(v,b,B,steps if i==0 else max(24,steps//96))
                    for i,(v,(b,B)) in enumerate(zip(durations,phases))]
        correction,_ = flow(initial,segments,sm["D"],sm["C"])
        return matrix.solve(arb_mat([[a+2*pi*c] for a,c in zip(base,correction)])).entries()
    out = image(L,t,2048)
    assert L.contains_interior(out[0]) and t.contains_interior(out[1]), (out[:2],L,t)
    for steps in (2048,4096,8192):
        L,t = L.intersection(out[0]),t.intersection(out[1])
        out = image(L,t,steps)
    L,t = L.intersection(out[0]),t.intersection(out[1])
    mu = model["mz"]*L.exp()
    kappa = mu*(-t).exp()
    # 最轻中心偶质量来自一个中性重矢量，倍率超过临界值即全线脱耦。
    even = [*layers[0],*(v for q in (4,2,6) for v in layers[q][3:])]
    floor = even[0]
    for v in even[1:]:
        floor = floor.min(v)
    assert floor.overlaps(layers[0][0]) and all(layers[q][3]>floor for q in (4,2,6))
    critical = t.exp()/floor
    assert kappa*min(v.lower() for q in (1,3) for v in layers[q])>model["mz"]
    assert 19 < critical < 23 and 38 < out[2] < 40
    assert arb("2e14")<mu<arb("5e14") and arb("2e12")<kappa<arb("6e12")
    print("十九标量" if full else "十六标量","激活的 Dirac 和标量数",
          sum(r[3] for r in active),sum(r[4] for r in active),
          "交点、质量单位/GeV、共同逆耦合、偶模脱耦临界倍率",
          *[(v.lower().str(9),v.upper().str(9)) for v in (mu,kappa,out[2],critical)])


def verify():
    model, _, _, initial, _ = chiral["source"]["mean"]["inputs"]()
    ctx.prec = 192
    read = chiral["source"]["spectral_supplier"](model)
    rows = read(acb.zeta_zero(11661).imag)
    layers = {q: chiral["ordered"]([v for v,z in rows if z==q]) for q in (0,*model["charges"])}
    assert all(layers[q][2]<9<layers[q][3] for q in model["charges"])
    for full in (True,False):
        crossing(model,layers,initial,full)


if __name__=="__main__":
    verify()
