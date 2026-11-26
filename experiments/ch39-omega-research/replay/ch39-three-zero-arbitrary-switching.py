"""三个非中心编号的作用在总时长不变时任意交错的检验。"""
from functools import lru_cache
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb

two=run_path(str(Path(__file__).with_name("ch39-two-zero-arbitrary-switching.py")))
three,free,ordered,study,box,support,edges=(two[k] for k in ("three","free","ordered","study","box","support","edges"))
flavor=ordered["flavor"]
whole=arb(0).union(1)


def prepare_component(values):
    lower=sorted(v.lower() for v in values)
    upper=sorted((v.upper() for v in values),reverse=True)
    def prefix(bounds):
        out=[arb(0)]
        for v in bounds:
            out.append(out[-1]+v/len(bounds))
        return out
    left,right=prefix(lower),prefix(upper)
    size=len(values)
    magnitude=max(abs(lower[0]),abs(upper[0]))
    def integral(duration):
        mid=duration.mid()
        count=min(size-1,int((mid*size).floor().unique_fmpz()))
        part=mid-arb(count)/size
        base=(left[count]+part*lower[count]).union(right[count]+part*upper[count])
        return base+arb(0,duration.rad()*magnitude)
    return integral


def verify():
    cost,_,scale,*_=support["prev"]["setup"]()
    cube=[v+arb(0,arb(2)**-28) for v in box["centers"]]
    basis=flavor["basis"]
    entries=free["entries"]
    radial=(1,5,29)
    attenuation=[((1-arb(radial[i])**2/867)*(1-arb(radial[j])**2/867)).sqrt() for i,j in entries]
    @lru_cache(None)
    def atom(n):
        height=acb.zeta_zero(n).imag
        x=arb(0) if n==11431 else study["center"]/(height**2+study["offset"])-1
        return x,arb(0) if n==11431 else cost(x)[0],study["native"]["plane_response"](height)
    @lru_cache(None)
    def reference(n):
        block=ordered["evolution"](atom(n)[2],arb(1),cube)
        wave=flavor["radial"]+basis*block*basis.transpose()
        return [wave[i,j] for i,j in entries],three["spectral_modes"](atom(n)[2],cube)
    @lru_cache(2048)
    def kernel(n,other,size):
        base,modes=reference(n)
        diff=cube[0]*(atom(other)[2]-atom(n)[2])
        norm=sum(abs(v)**2 for v in diff.entries()).sqrt()
        terms=[]
        for a,qa in modes:
            for b,qb in modes:
                mat=acb(0,1)*acb(0,a).exp()*basis*qa*diff*qb*basis.transpose()
                terms.append((mat,b-a))
        columns=[[] for _ in entries]
        for k in range(size):
            t=arb(k).union(k+1)/size
            val=acb_mat(3,3)
            for matrix,speed in terms:
                val+=acb(0,speed*t).exp()*matrix
            for accum,(i,j) in zip(columns,entries):
                accum.append(val[i,j])
        prepared=[(prepare_component([v.real for v in col]),prepare_component([v.imag for v in col]))
                  for col in columns]
        def integral(duration):
            return [acb(re(duration),im(duration)) for re,im in prepared]
        return norm,integral
    def ruled(n,shorts,durations,extra=0,size=8):
        base,_=reference(n)
        values=list(base)
        norm=arb(0)
        for k,p in zip(shorts,durations):
            p=p.intersection(whole)
            bound,integral=kernel(n,k,size)
            values=[a+b for a,b in zip(values,integral(p))]
            norm+=p*bound
        error=norm**2/2+extra
        padded=[v+acb(arb(0,e*error),arb(0,e*error)) for v,e in zip(values,attenuation)]
        return free["excluded"](padded)
    corner=[v+s*arb(2)**-28 for v,s in zip(box["centers"],(1,-1,1,-1))]
    data=list(map(atom,(11431,11432,11415,11430,11372,11304)))
    weights,_=edges["candidate"]([x+scale*c/2 for x,c,t in data],
                                 [box["flavor_data"](t,corner)[0] for x,c,t in data],arb(1)/3)
    ceiling=support["dot"](weights,[c for x,c,t in data])
    eps=scale*ceiling/2
    gauge=(-eps).union(eps)
    max_contrast=abs(cube[0])*2/(arb("1.3")-4*arb.pi()/32)
    def lower(l,h):
        x,c,_=atom(l)
        y,d,_=atom(h)
        return d+(gauge-y)/(x-y)*(c-d)
    assert atom(11430)[0]>eps and atom(11432)[0]<-eps
    assert lower(11390,11432)>ceiling and lower(11430,11500)>ceiling
    def examine(l,h):
        x,c,_=atom(l)
        y,d,_=atom(h)
        p=(gauge-y)/(x-y)
        extra=(ceiling-lower(l,h)).upper()
        back,short=(h,l) if p<arb(1)/2 else (l,h)
        q=p if short==l else 1-p
        def parameters(n):
            v,cn,_=atom(n)
            k=(v-y)/(x-y)
            residual=cn-(d+k*(c-d))
            assert residual>0
            return k,(extra/residual).upper()
        def tail(n):
            k,tmax=parameters(n)
            # 合入背景作用以消去远处编号；所有短态时长改变由弦斜率界控制。
            spread=abs(k if short==l else 1-k)*tmax
            return ruled(back,[short],[q+arb(0,spread)],extra=tmax*max_contrast)
        cuts=[]
        for sign,base in ((-1,l),(1,h)):
            step=8
            while not tail(base+sign*step):
                step*=2
                assert step<=8192
            cuts.append(base+sign*step)
        cells=0
        for n in (*range(cuts[0]+1,l),*range(h+1,cuts[1])):
            k,cap=parameters(n)
            last=min(cap,(p/k if n<l else (1-p)/(1-k)).upper())
            queue=[(arb(0),last,0)]
            while queue:
                lo,hi,depth=queue.pop()
                t=lo.union(hi)
                durations=[p-k*t,1-p+(k-1)*t,t]
                ns=[l,h,n]
                background=max(range(3),key=lambda j:float(durations[j].mid()))
                chosen=ns.pop(background)
                durations.pop(background)
                cells+=1
                if ruled(chosen,ns,durations,size=8):
                    continue
                if depth>=4 and ruled(chosen,ns,durations,size=64):
                    continue
                assert depth<14,(l,h,n,t.str(14))
                mid=(lo+hi)/2
                queue.extend(((lo,mid,depth+1),(mid,hi,depth+1)))
        print("Pair",l,h,"tails",cuts,"occupation cells",cells,flush=True)
        return cells
    totals=[]
    for l in range(11391,11431):
        for h in range(11432,11500):
            value=lower(l,h)
            if value>ceiling:
                break
            assert value<ceiling
            totals.append(examine(l,h))
    assert len(totals)==20
    print("All center-free controls on up to three generators excluded",sum(totals))


if __name__=="__main__":
    verify()
