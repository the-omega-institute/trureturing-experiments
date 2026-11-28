"""源 2345698 的四条双支点线逐一独立调整，认证引力根的全域唯一性及规范残差。"""
from math import factorial
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

ROOT = Path(__file__).resolve().parent
TARGETS = (*range(2345700,2345747,2),*range(2345964,2345987,2))


def verify():
    repair = run_path(str(ROOT/"ch39-one-share-gravity-repair.py"))
    load = repair["load"]
    edges,kernel,moments,lines = repair["source_data"]()
    _,lib,spectra = load("two-zero-cutoff-gravity")["gravity_fields"](
        ROOT,(repair["SOURCE"],*TARGETS),tight_neutral=True)
    model,_,_,initial,_ = lib["chiral"]["source"]["mean"]["inputs"]()
    fields = [(0,0,0),*load("signed-gravity-unification")["gauge_weights"](model)]
    beta = model["beta"]
    t = (beta[0]-beta[2])/(beta[0]-beta[1])
    q = [1-t,t,-1]
    baseline = sum(c*a for c,a in zip(q,initial,strict=True))
    norm = abs(q[0])+abs(q[1])
    gauge = arb_mat([[0,0,*q]])
    scale = arb(7)/1000
    weights = sum(abs(w) for a,b,w in edges)
    order = kernel["ORDER"]+1
    factor = 2/(4*arb.pi()*factorial(order))
    min_x = min(a.lower() for a,b,w in edges)*scale/10
    assert all(max(z.upper() for g,a,z in rows)<min(arb(1),2*min_x) for rows in spectra)
    envelope = [(g,a,max(bg[i][2].upper() for bg in spectra[1:])) for i,(g,a,z) in enumerate(spectra[0])]
    def loss(rows,s):
        def at(b):
            x = scale*s*b
            return x*sum(abs(a)*(1-lib["expint"](2,z/x)) for g,a,z in rows)
        return sum(abs(w)*at(b) for _,b,w in edges)/(2*arb.pi())
    source_loss,target_loss = loss(spectra[0],1),loss(envelope,arb(1)/10)
    massless = sum(a for g,a,z in spectra[0])/(2*arb.pi())
    divide = load("one-share-gravity-unique")["divided"]
    source_cut = arb(".99")*scale*min(a.lower() for a,b,w in edges)
    identity = arb_mat([[1,0],[0,1]])
    gauge_tail = factor*weights*sum(sum(abs(c*b) for c,b in zip(q,bs,strict=True)) for bs in fields)
    global_gaps = []
    for number,rows in zip(TARGETS,spectra[1:],strict=True):
        mat = kernel["series"](rows,fields)
        coeff = [arb_mat([[int(i==k) for i in range(5)]])*mat for k in (0,1)]
        tails = [factor*sum(abs(2*a) for g,a,z in rows),factor*sum(abs(g) for g,a,z in rows)]
        assert max(z.upper() for g,a,z in rows)<source_cut
        gaps = []
        for line,(change,share) in lines.items():
            physical = (-share).union(1-share)
            small = (1+source_loss+abs(physical*massless*change[1,0]*kernel["S"]))/10+target_loss
            assert small<1,(line,small)
            args = (coeff,moments,change,[t*weights for t in tails],[2*t*abs(edges[line][2]) for t in tails])
            root,_ = repair["repair_box"](args,share)
            s0,d0 = root
            # 在唯一性盒内对四种不同大小的斜率分别做预条件。
            radius = arb("0.00002")
            domain = [center+arb(0,radius) for center in root]
            assert domain[0]>arb(".99") and 0<share+domain[1]<1
            _,jac = repair["bounds"](*domain,*args)
            _,center = repair["bounds"](*root,*args)
            inv = arb_mat(2,2,[v.mid() for v in center.inv().entries()])
            assert not inv.det().contains(0) and max(repair["slopes"](identity-inv*jac))<1
            left = arb(1)/10
            pending = []
            while left<1/scale:
                pending.append((left,(2*left).min(1/scale)))
                left *= 2
            while pending:
                lo,hi = pending.pop()
                box = lo.union(hi)
                kn,wn = [divide(row,box,i+1) for i,row in enumerate(coeff)]
                slope = (wn*change)[0,0]
                assert not slope.contains(0)
                werr = (tails[1]*weights*2**order/box**2).upper()
                d = (arb("1e-122")/box**2-(wn*moments)[0,0]+arb(0,werr))/slope
                if not physical.overlaps(d):
                    continue
                d = physical.intersection(d)
                kerr = (tails[0]*weights*2**order/box).upper()
                inv_s = (kn*(moments+d*change))[0,0]+arb(0,kerr)
                if not inv_s.overlaps(1/box):
                    continue
                possible = 1/inv_s.intersection(1/box)
                if domain[0].contains(possible) and domain[1].contains(d):
                    continue
                mid = ((lo+hi)/2).mid()
                assert hi-lo>arb("1e-12"),(number,line,domain,d,possible)
                pending.extend(((lo,mid),(mid,hi)))
            score = baseline+(repair["dilate"](gauge*mat,s0)*(moments+d0*change))[0,0]+arb(0,gauge_tail.upper())
            gap = (abs(score)/norm).lower()
            assert gap>0
            gaps.append(gap)
        known = {2345704:".03667",2345712:".01346",2345716:".04801",2345972:".04281",2345974:".04803"}
        if number in known:
            assert min(gaps)>arb(known[number])
        global_gaps.append(min(gaps))
        print(number,"四条线按 16、285、301、323 的最大绝对规范差下界",*[g.str(12) for g in gaps],flush=True)

    assert len(global_gaps)==36 and min(global_gaps)>arb(".01346")
    print("全部单线组合的规范差统一下界",min(global_gaps).str(14),flush=True)


if __name__=="__main__":
    verify()
