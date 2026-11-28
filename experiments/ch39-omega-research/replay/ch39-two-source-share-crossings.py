"""两条原有谱线只换各自两个原支点间的份额，认证五个重谱的正交汇。"""
import json
from math import factorial
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

ROOT = Path(__file__).resolve().parent


def verify():
    repair = run_path(str(ROOT/"ch39-one-share-gravity-repair.py"))
    load = repair["load"]
    edges,kernel,moments,choices = repair["source_data"]()
    seeds = json.loads((ROOT/"ch39-two-source-share-seeds.json").read_text())
    targets = sorted({r["heavy"] for r in seeds})
    assert targets == [2345704,2345712,2345716,2345972,2345974]
    assert len(seeds)==len({(r["heavy"],*r["lines"]) for r in seeds})==8
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
            def read(s,d,e):
                measure = moments+d*mi+e*mj
                residual, jac = [], []
                for k,bias in enumerate((-1,-arb("1e-122"),constant)):
                    r = row(coeff,k)
                    f = repair["dilate"](r,s)
                    derivative = repair["dilate"](r,s,True)
                    tail = factor*size[k]
                    residual.append((f*measure)[0,0]+bias+arb(0,(tail*total).upper()))
                    jac.append([(derivative*measure)[0,0]+arb(0,(2*tail*total).upper()),
                                (f*mi)[0,0]+arb(0,(2*tail*abs(edges[i][2])).upper()),
                                (f*mj)[0,0]+arb(0,(2*tail*abs(edges[j][2])).upper())])
                return arb_mat([[v] for v in residual]),arb_mat(jac)
            center = [arb(v)/10**8 for v in seed["center"]]
            for _ in range(2):
                value,jac = read(*center)
                step = jac.solve(value)
                center = [(x-y).mid() for x,y in zip(center,step.entries(),strict=True)]
            radius = arb("1e-9")
            box = [x+arb(0,radius) for x in center]
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
            measure = moments+d*mi+e*mj
            coupling = []
            for k in range(3):
                f = repair["dilate"](row(series,2+k),s)
                tail = factor*total*sum(abs(bs[k]) for bs in fields)
                coupling.append(initial[k]+(f*measure)[0,0]+arb(0,tail.upper()))
            log_energy = (coupling[0]-coupling[2])/(beta[0]-beta[2])
            energy = model["mz"]*log_energy.exp()
            common = coupling[2]-beta[2]*log_energy
            assert model["mz"]<energy<arb("1.22089e19") and min(*coupling,common)>arb("19.83")
            if number==2345704 and [i,j]==[16,301]:
                assert arb("158493.00") < energy < arb("158493.02")
            print(number,i,j,"标尺比、两份额增量、交汇 GeV、全程最小逆耦合",
                  *[v.str(9) for v in (s,d,e,energy,min(*coupling,common))],flush=True)


if __name__=="__main__":
    verify()
