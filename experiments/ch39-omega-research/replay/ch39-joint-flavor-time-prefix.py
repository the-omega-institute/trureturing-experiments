"""三代、两圈与曲率窗口共用参考味尺度时的算术前缀。"""
from concurrent.futures import ProcessPoolExecutor
from functools import cache
from pathlib import Path
from runpy import run_path

from flint import acb,arb,ctx


root=Path(__file__).parent
source=run_path(str(root/"ch39-flavor-scale-time-prefix.py"))
physical=run_path(str(root/"ch39-first-three-mode-observables.py"))
shell=run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
indices=(2529893,2562426,4035321,4193973,6736068,6964994,25215595,25408663)


@cache
def physical_tests():
    ctx.prec=112
    accepted=physical["mask"]["three_mode_test"]()
    coefficient=(9*arb(3).log()-12*arb(2).log())/2
    light_scale=arb("246.22")**2/(32*arb("2.435e18")/724**2*arb("1.22089e19"))
    return accepted,coefficient,light_scale


def check_gap_block(block):
    first,count=block
    test,coefficient,light_scale=physical_tests()
    ctx.prec=80
    zeros=acb.zeta_zeros(first,count)
    assert len(zeros)==count
    three_modes=full_solves=joint_hits=0
    for n,z in enumerate(zeros,first):
        ctx.prec=112
        h=z.imag
        if not test(h):
            continue
        three_modes+=1
        delta=arb(23)/25*arb(2)**-20/(h**2+shell["soft"]["point"]["study"]["offset"])
        center=15-delta/3
        eig=shell["laplacian"](2,h).eig(multiple=True)
        assert len(eig)==19 and all(v.imag.contains(0) for v in eig)
        mass=physical["joint"]["chiral"]["ordered"]([abs(v.real*(v.real-center)) for v in eig])[15]
        a=(light_scale/mass)**2
        curve=a**2*coefficient/(arb.pi()*arb("1e-122"))*(1-arb(0).union(a/(arb(".03")*coefficient)))
        if abs(curve-1)>arb(".002"):
            continue
        full_solves+=1
        _,alpha,curve=physical["observables"](n)
        if physical["inside"](alpha,arb(".117"),arb(".119")) and physical["inside"](curve,arb(".998"),arb("1.002")):
            joint_hits+=1
    return count,three_modes,full_solves,joint_hits


def verify():
    left,right=indices[:2]
    blocks=[(n,min(768,right-n)) for n in range(left+1,right,768)]
    with ProcessPoolExecutor(max_workers=4) as workers:
        totals=tuple(map(sum,zip(*workers.map(check_gap_block,blocks),strict=True)))
        responses=list(workers.map(physical_worker,indices))
    assert totals[0]==right-left-1 and totals[3]==0
    print("间隙内枚数、三模数、完整交点复算数、共同命中数",totals)
    ctx.prec,ctx.cap=1536,3
    times=source["clock"]["clock"]["times"]()
    lo,hi=source["determinant"]["scale_interval"]()
    scale=lo.union(hi)
    ctx.prec=256
    heights=[acb.zeta_zero(n).imag for n in indices]
    assert heights[0]<scale*acb.zeta_zero(1).imag<heights[1]
    ctx.prec=1536
    weights=source["certify_weights"](times,scale,7,heights)
    margin=source["separation_margin"](times,scale,8,0,heights[0],heights[1])
    assert margin>arb("2.43e-18")
    assert min(p.lower() for p in weights)>arb(".0213")
    probability=sum(p*arb(vub)**2 for p,vub in zip(weights,responses,strict=True))
    assert arb("8.86770e-6")<probability<arb("8.86772e-6")
    print("八项最大误差下界、未分辨编号的混合概率",margin.str(12),probability.str(12))


def physical_worker(n):
    _,alpha,curve=physical["observables"](n)
    assert arb(".117")<alpha<arb(".119") and arb(".998")<curve<arb("1.002")
    flavor=physical["joint"]["flavor"]
    vub=flavor["reference_readout"]()(flavor["singlet"](n)[0])[2]
    return vub.str(24,more=True)


if __name__=="__main__":
    verify()
