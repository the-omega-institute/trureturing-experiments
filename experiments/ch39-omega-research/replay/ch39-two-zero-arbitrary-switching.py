"""检验只用两个非中心生成元任意切换的规范、端口和 CKM 条件。"""
from functools import lru_cache
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb

root=Path(__file__).parent
free=run_path(str(root/"ch39-center-free-three-pulses.py"))
three,ordered,study,box,support,edges=(free[k] for k in ("three","ordered","study","box","support","edges"))
flavor=ordered["flavor"]


def interval_integral(bounds,duration):
    # 任意可测的作用时段，总时长已知。各格至多占满一格。
    size=len(bounds)
    count=int((size*duration).floor().unique_fmpz())
    remainder=duration-arb(count)/size
    mins=sorted(v.lower() for v in bounds)
    maxs=sorted((v.upper() for v in bounds),reverse=True)
    low=sum(mins[:count])/size+remainder*mins[count]
    high=sum(maxs[:count])/size+remainder*maxs[count]
    return low.union(high)


def check_switches(short,long,duration,coeff,cells):
    basis=flavor["basis"]
    base=ordered["evolution"](long,arb(1),coeff)
    reference=flavor["radial"]+basis*base*basis.transpose()
    diff=coeff[0]*(short-long)
    # 两个本征投影将相互作用表象中的一阶积分写成常数加两个有界相位项。
    modes=three["spectral_modes"](long,coeff)
    gaps=[]
    for a,qa in modes:
        for b,qb in modes:
            matrix=acb(0,1)*basis*qa*diff*qb*basis.transpose()
            gaps.append((acb(0,a).exp()*matrix,b-a))
    integrands=[[] for _ in free["entries"]]
    for i in range(cells):
        t=arb(i).union(i+1)/cells
        change=acb_mat(3,3)
        for matrix,gap in gaps:
            change+=acb(0,t*gap).exp()*matrix
        for col,(k,j) in zip(integrands,free["entries"]):
            col.append(change[k,j])
    radial=(1,5,29)
    delta=duration*sum(abs(x)**2 for x in diff.entries()).sqrt()
    remainder=delta**2/2
    values=[]
    for col,(k,j) in zip(integrands,free["entries"]):
        linear=acb(interval_integral([v.real for v in col],duration),
                   interval_integral([v.imag for v in col],duration))
        error=remainder*((1-arb(radial[k])**2/867)*(1-arb(radial[j])**2/867)).sqrt()
        values.append(reference[k,j]+linear+acb(arb(0,error),arb(0,error)))
    return free["excluded"](values)


def verify():
    cost,_,scale,*_=support["prev"]["setup"]()
    cube=[v+arb(0,arb(2)**-28) for v in box["centers"]]
    @lru_cache(None)
    def atom(n):
        height=acb.zeta_zero(n).imag
        x=arb(0) if n==11431 else study["center"]/(height**2+study["offset"])-1
        return x,arb(0) if n==11431 else cost(x)[0],study["native"]["plane_response"](height)
    corner=[v+s*arb(2)**-28 for v,s in zip(box["centers"],(1,-1,1,-1))]
    atoms=[atom(n) for n in (11431,11432,11415,11430,11372,11304)]
    w,_=edges["candidate"]([x+scale*c/2 for x,c,t in atoms],
                           [box["flavor_data"](t,corner)[0] for x,c,t in atoms],arb(1)/3)
    ceiling=support["dot"](w,[c for x,c,t in atoms])
    eps=scale*ceiling/2
    gauge=(-eps).union(eps)
    def lower(l,h):
        x,c,_=atom(l)
        y,d,_=atom(h)
        return d+(gauge-y)/(x-y)*(c-d)
    assert atom(11430)[0]>eps and atom(11432)[0]<-eps
    assert lower(11390,11432)>ceiling and lower(11430,11500)>ceiling
    attempts=[]
    for l in range(11391,11431):
        for h in range(11432,11500):
            level=lower(l,h)
            if level>ceiling:
                break
            assert level<ceiling
            x,c,left=atom(l)
            y,d,right=atom(h)
            duration=(gauge-y)/(x-y)
            assert 0<duration<1
            if duration>arb(1)/2:
                duration=1-duration
                left,right=right,left
            for cells in (8,32,128,512):
                if check_switches(left,right,duration,cube,cells):
                    attempts.append(cells)
                    break
            else:
                print("Unresolved",l,h,flush=True)
    print("Certified pairs",len(attempts),"cell counts",attempts)
    assert len(attempts)==20


if __name__=="__main__":
    verify()
