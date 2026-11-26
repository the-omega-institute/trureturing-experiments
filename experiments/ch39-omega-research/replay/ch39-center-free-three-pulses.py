"""检验不含 11431 的三个有序零点作用是否能进入联合误差预算。"""
from functools import lru_cache
from itertools import permutations, product
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

root=Path(__file__).parent
budget=run_path(str(root/"ch39-global-three-zero-clock.py"))
three=run_path(str(root/"ch39-ordered-three-zero-minimum.py"))
ordered,study,box,support,edges=(budget[k] for k in ("ordered","study","box","support","edges"))
flavor=ordered["flavor"]
entries=((0,1),(1,2),(0,2),(1,1))
tolerances=[arb(".000055")]*3+[arb(".006")]


def excluded(values):
    us,cb,ub,cs=values
    a,b,u=[abs(z)**2 for z in (us,cb,ub)]
    if u>arb(".003787")**2 or u<arb(".003677")**2:
        return True
    if not (1-u)>0:
        return False
    if any((prob-arb(lo)**2*(1-u)<0 or prob-arb(hi)**2*(1-u)>0)
           for prob,lo,hi in ((a,".224955",".225065"),(b,".041775",".041885"))):
        return True
    phase=us*cb*ub.conjugate()*cs.conjugate()+a*b*u/(1-u)
    if phase.real<0 or phase.imag<0:
        return True
    return phase.imag-arb("1.141").tan()*phase.real<0 or phase.imag-arb("1.153").tan()*phase.real>0


def verify():
    cost,_,scale,*_=support["prev"]["setup"]()
    theta=arb(1)/2
    cube=[c+arb(0,arb(2)**-28) for c in box["centers"]]
    @lru_cache(None)
    def coordinates(n):
        height=acb.zeta_zero(n).imag
        x=arb(0) if n==11431 else study["center"]/(height**2+study["offset"])-1
        return height,x,arb(0) if n==11431 else cost(x)[0]
    @lru_cache(None)
    def response(n):
        return study["native"]["plane_response"](coordinates(n)[0])
    @lru_cache(None)
    def modes(n):
        return three["spectral_modes"](response(n),cube)
    corner=[v+s*arb(2)**-28 for v,s in zip(box["centers"],(1,-1,1,-1))]
    ids=(11431,11432,11415,11430,11372,11304)
    vals=[coordinates(n) for n in ids]
    ws,_=edges["candidate"]([x+theta*scale*c for h,x,c in vals],
                            [box["flavor_data"](response(n),corner)[0] for n in ids],arb(1)/3)
    ceiling=support["dot"](ws,[c for h,x,c in vals])
    eps=scale*ceiling/2
    # 已证的全相位内块下界 1.3-4π/32 与两列投影的算子范数平方上界 2。
    # 两个响应都半负定，所以算子差的范数受同一个正上界控制。
    max_contrast=abs(cube[0])*2/(arb("1.3")-4*arb.pi()/32)
    radial=(1,5,29)
    attenuation=[((1-arb(radial[i])**2/867)*(1-arb(radial[j])**2/867)).sqrt() for i,j in entries]
    gauge=(-eps).union(eps)
    basis=flavor["basis"]
    def paired(l,h):
        _,x,c=coordinates(l)
        _,y,d=coordinates(h)
        mass=(gauge-y)/(x-y)
        return x,y,c,d,mass
    def min_cost(l,h):
        x,y,c,d,p=paired(l,h)
        return d+p*(c-d)
    assert coordinates(11430)[1]>eps and coordinates(11432)[1]<-eps
    assert min_cost(11390,11432)>ceiling and min_cost(11430,11500)>ceiling
    def study_pair(l,h):
        x,y,c,d,p=paired(l,h)
        extra=(ceiling-(d+p*(c-d))).upper()
        assert extra>0
        tl,th=response(l),response(h)
        diff=tl-th
        contrast=abs(cube[0])*sum(abs(v)**2 for v in diff.entries()).sqrt()
        pair_projection=[]
        for order in ((l,h),(h,l)):
            steps={l:ordered["evolution"](tl,p,cube),h:ordered["evolution"](th,1-p,cube)}
            reduced=steps[order[0]]*steps[order[1]]
            ref=flavor["radial"]+basis*reduced*basis.transpose()
            pair_projection.append([ref[i,j] for i,j in entries])
        def near_pair(error):
            padding=[acb(arb(0,error*s),arb(0,error*s)) for s in attenuation]
            return all(excluded([v+e for v,e in zip(base,padding)]) for base in pair_projection)
        def tail(n):
            _,v,b=coordinates(n)
            k=(v-y)/(x-y)
            remainder=b-(d+k*(c-d))
            assert remainder>0
            mass=extra/remainder
            error=(max_contrast+contrast*max(abs(k),abs(1-k)))*mass
            return near_pair(error)
        cuts=[]
        for sign,base in ((-1,l),(1,h)):
            step=8
            while not tail(base+sign*step):
                step*=2
                assert step<=8192,(l,h,sign)
            cuts.append(base+sign*step)
        cells=0
        triples=0
        for n in (*range(cuts[0]+1,l+1),*range(h,cuts[1])):
            # 包含重复编号，此时时长变量是在同一编号的两次作用之间分配。
            if n==l:
                k=arb(1)
                last=p.upper()
            elif n==h:
                k=arb(0)
                last=(1-p).upper()
            else:
                _,v,b=coordinates(n)
                k=(v-y)/(x-y)
                residual=b-(d+k*(c-d))
                assert residual>0
                by_positive=p/k if n<l else (1-p)/(1-k)
                last=min((extra/residual).upper(),by_positive.upper())
            assert last>0
            tk=response(n)
            local=max(sum(abs(v)**2 for v in (tk-t).entries()).sqrt() for t in (tl,th))
            error=(abs(cube[0])*local+contrast*max(abs(k),abs(1-k)))*last
            if near_pair(error):
                continue
            triples+=1
            ns=(l,h,n)
            for order in permutations(range(3)):
                terms=[]
                for signs in product((0,1),repeat=3):
                    proj=flavor["identity"]
                    for i in order:
                        proj=proj*modes(ns[i])[signs[i]][1]
                    coeff=basis*proj*basis.transpose()
                    lam=[modes(label)[sgn][0] for label,sgn in zip(ns,signs)]
                    start=lam[1]+p*(lam[0]-lam[1])
                    speed=(lam[2]-lam[1])-k*(lam[0]-lam[1])
                    terms.append((start,speed,[coeff[i,j] for i,j in entries]))
                slopes=[sum(abs(speed)*abs(col[j]) for _,speed,col in terms) for j in range(4)]
                queue=[(arb(0),last)]
                while queue:
                    lo,hi=queue.pop()
                    mid=(lo+hi)/2
                    width=(hi-lo)/2
                    values=[flavor["radial"][i,j] for i,j in entries]
                    for start,speed,col in terms:
                        phase=acb(0,start+speed*mid).exp()
                        values=[v+phase*c for v,c in zip(values,col)]
                    values=[v+acb(arb(0,width*s),arb(0,width*s)) for v,s in zip(values,slopes)]
                    cells+=1
                    if excluded(values):
                        continue
                    assert width>arb("1e-10"),(ns,order,mid.str(18))
                    queue.extend(((lo,mid),(mid,hi)))
        print("Pair",l,h,"tail boundaries",*cuts,"explicit triples",triples,"time cells",cells,flush=True)
        return triples,cells
    tested=[]
    for l in range(11391,11431):
        for h in range(11432,11500):
            low=min_cost(l,h)
            if low>ceiling:
                break
            assert low<ceiling
            tested.append(study_pair(l,h))
    print("Center-free ordered solutions excluded",len(tested),*[sum(v) for v in zip(*tested)],flush=True)


if __name__=="__main__":
    verify()
