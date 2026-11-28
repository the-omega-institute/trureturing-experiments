"""两条原有谱线只换各自两个原支点间的份额，认证五个重谱的正交汇。"""
import json
from math import factorial
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

ROOT = Path(__file__).resolve().parent


FORWARD = tuple(range(2_345_988, 2_346_027, 2))


def certified_pairs(targets,seeds):
    """逐个给出三元根、全程读数及区间残差函数。"""
    repair = run_path(str(ROOT/"ch39-one-share-gravity-repair.py"))
    load = repair["load"]
    edges,kernel,moments,choices = repair["source_data"]()
    _,lib,spectra = load("two-zero-cutoff-gravity")["gravity_fields"](ROOT,targets,tight_neutral=True)
    model,_,_,initial,_ = lib["chiral"]["source"]["mean"]["inputs"]()
    fields = [(0,0,0),*load("signed-gravity-unification")["gauge_weights"](model)]
    beta = [b/(2*arb.pi()) for b in model["beta"]]
    t = (beta[0]-beta[2])/(beta[0]-beta[1])
    direction = [1-t,t,-1]
    constant = sum(v*a for v,a in zip(direction,initial,strict=True))
    matrix = arb_mat([[1,0,0,0,0],[0,1,0,0,0],[0,0,*direction]])
    total = sum(abs(w) for _,_,w in edges)
    factor = 2/(4*arb.pi()*factorial(kernel["ORDER"]+1))
    def row(mat,i):
        return arb_mat([[mat[i,j] for j in range(mat.ncols())]])
    for number,rows in zip(targets,spectra,strict=True):
        assert max(z.upper() for _,_,z in rows)<min(arb(1),arb(".0069")*min(a.lower() for a,_,_ in edges))
        series = kernel["series"](rows,fields)
        coeff = matrix*series
        size = [sum(abs(2*a) for _,a,_ in rows),sum(abs(g) for g,_,_ in rows),
                sum(abs(d*b) for bs in fields for d,b in zip(direction,bs,strict=True))]
        for seed in (r for r in seeds if r["heavy"]==number):
            i,j = seed["lines"]
            mi,pi = choices[i]
            mj,pj = choices[j]
            scale = arb(7)/1000
            def variation(d,e):
                return total+2*(abs(d*edges[i][2])+abs(e*edges[j][2]))
            def observable(s,d,e):
                assert max(z.upper() for _,_,z in rows) < s*scale*min(a.lower() for a,_,_ in edges)
                measure = moments+d*mi+e*mj
                coupling = []
                for k in range(3):
                    f = repair["dilate"](row(series,2+k),s)
                    tail = factor*variation(d,e)*sum(abs(bs[k]) for bs in fields)
                    coupling.append(initial[k]+(f*measure)[0,0]+arb(0,tail.upper()))
                log_energy = (coupling[0]-coupling[2])/(beta[0]-beta[2])
                energy = model["mz"]*log_energy.exp()
                common = coupling[2]-beta[2]*log_energy
                positions = (pi+d,pj+e)
                floor = min(*coupling,common)
                checks = [s,1-scale*s,log_energy,(arb("1.22089e19")/model["mz"]).log()-log_energy,
                          *positions,*[1-v for v in positions],floor]
                assert all(not v.contains(0) for v in checks)
                labels = ("正标尺","标尺上限","能标下限","能标上限","第一份额下限","第二份额下限",
                          "第一份额上限","第二份额上限","正性")
                failed = [label for label,v in zip(labels,checks,strict=True) if v<0]
                return energy,floor,failed
            def read(s,d,e):
                assert s>arb(".5")
                assert max(z.upper() for _,_,z in rows) < s*scale*min(a.lower() for a,_,_ in edges)
                measure = moments+d*mi+e*mj
                residual, jac = [], []
                for k,bias in enumerate((-1,-arb("1e-122"),constant)):
                    r = row(coeff,k)
                    f = repair["dilate"](r,s)
                    derivative = repair["dilate"](r,s,True)
                    tail = factor*size[k]
                    residual.append((f*measure)[0,0]+bias+arb(0,(tail*variation(d,e)).upper()))
                    jac.append([(derivative*measure)[0,0]+arb(0,(2*tail*variation(d,e)).upper()),
                                (f*mi)[0,0]+arb(0,(2*tail*abs(edges[i][2])).upper()),
                                (f*mj)[0,0]+arb(0,(2*tail*abs(edges[j][2])).upper())])
                return arb_mat([[v] for v in residual]),arb_mat(jac)
            approximate = "center" not in seed
            proposal = penergy = predicted = None
            if approximate:
                source = [arb(1),arb(0),arb(0)]
                residual,jac = read(*source)
                delta = jac.solve(residual)
                proposal = [x-y for x,y in zip(source,delta.entries(),strict=True)]
                penergy,_,pfail = observable(*proposal)
                predicted = not pfail
                center = [x.mid() for x in proposal]
            else:
                center = [arb(v)/10**8 for v in seed["center"]]
            for _ in range(4 if approximate else 2):
                value,jac = read(*center)
                step = jac.solve(value)
                center = [(x-y).mid() for x,y in zip(center,step.entries(),strict=True)]
            radius = arb("1e-9")
            box = [x+arb(0,radius) for x in center]
            if not approximate:
                assert box[0]>arb(".99") and box[0]<arb("1.02")
                assert 0<pi+box[1]<1 and 0<pj+box[2]<1
            value,jac = read(*center)
            inv = arb_mat(3,3,[x.mid() for x in jac.inv().entries()])
            assert not inv.det().contains(0)
            _,span = read(*box)
            resid = arb_mat([[int(a==b) for b in range(3)] for a in range(3)])-inv*span
            contraction = max(sum(abs(resid[a,b]) for b in range(3)).upper() for a in range(3))
            error = max(abs(z).upper() for z in (inv*value).entries())
            assert contraction<1 and error+radius*contraction < radius
            s,d,e = [x+arb(0,(error/(1-contraction)).upper()) for x in center]
            energy,minimum,failed = observable(s,d,e)
            yield number,(i,j),(s,d,e),energy,minimum,failed,read,proposal,penergy,predicted


def verify(*, forward=False, check_domain=None, targets_override=None):
    if targets_override is not None and not forward:
        raise ValueError("自选目标适用于从源点预报的检验")
    seeds = json.loads((ROOT/"ch39-two-source-share-seeds.json").read_text())
    targets = sorted({r["heavy"] for r in seeds})
    assert targets == [2345704,2345712,2345716,2345972,2345974]
    assert len(seeds)==len({(r["heavy"],*r["lines"]) for r in seeds})==8
    if forward:
        targets = list(FORWARD if targets_override is None else targets_override)
        seeds = [{"heavy":n, "lines":[16,301]} for n in targets]
    outcomes = []
    certified = {}
    for number,(i,j),root,energy,minimum,failed,read,proposal,penergy,predicted in certified_pairs(targets,seeds):
        s,d,e = root
        physical = not failed
        if check_domain is not None:
            check_domain(number,(s,d,e),physical,read)
        if forward:
            print(number,"一阶预报、认证分支、预报与实际能标、份额及标尺的最大坐标误差",
                  predicted,physical,penergy.str(9),energy.str(9),",".join(failed),
                  max(abs(x-y).upper() for x,y in zip(proposal,(s,d,e),strict=True)).str(9),flush=True)
            certified[number] = (physical,energy)
            if physical and targets_override is None:
                assert abs(penergy/energy-1) < arb(".00060")
                assert max(abs(x-y).upper() for x,y in zip(proposal,(s,d,e),strict=True)) < arb(".000055")
            outcomes.append((number,predicted,physical,abs((penergy/energy).log()).upper()))
            continue
        assert physical and minimum>arb("19.83")
        if number==2345704 and [i,j]==[16,301]:
            assert arb("158493.00") < energy < arb("158493.02")
        print(number,i,j,"标尺比、两份额增量、交汇 GeV、全程最小逆耦合",
              *[v.str(9) for v in (s,d,e,energy,minimum)],flush=True)
    if targets_override is not None:
        return certified
    if forward:
        assert [n for n,p,v,_ in outcomes if p and v] == [2345990,2345998,2346002,2346024,2346026]
        assert all(p==v for _,p,v,_ in outcomes) and max(err for _,_,_,err in outcomes)<arb(".06")
        print("二十处算术线性预报和数学根均正交汇",
              [n for n,p,v,_ in outcomes if p and v],flush=True)
        print("判定一致数，最大能标对数误差",sum(p==v for _,p,v,_ in outcomes),
              max(err for _,_,_,err in outcomes).str(9),flush=True)


if __name__=="__main__":
    verify()
    verify(forward=True)
