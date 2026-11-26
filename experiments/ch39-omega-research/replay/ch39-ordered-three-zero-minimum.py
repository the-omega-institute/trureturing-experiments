"""中心占三分之二且规范均值精确时，三段有序演化的全谱端口最小值。"""
from itertools import permutations, product
from pathlib import Path
from runpy import run_path

from flint import acb, arb

ordered = run_path(str(Path(__file__).with_name("ch39-ordered-flavor-unitary.py")))
study,flavor = ordered["study"],ordered["flavor"]
center = 11431


def weights(x,y,remaining):
    left = -remaining*y/(x-y)
    return (1-remaining,left,remaining-left)


def spectral_modes(response,coefficients):
    eta,h32,h152,alpha=coefficients
    ident=flavor["identity"]
    action=eta*response+h32*flavor["arcs"][32]+h152*flavor["arcs"][152]+alpha*ident
    mean=action.trace().real/2
    traceless=action-mean*ident
    gap=(traceless[0,0].real**2+abs(traceless[0,1])**2).sqrt()
    return [(mean+sgn*gap,(ident+sgn*traceless/gap)/2) for sgn in (-1,1)]


def exclusion(order,atoms,modes,shares,end,tolerances,start=None):
    # 当前外部质量为 r，各段时长为 (1-r, r w, r(1-w))。
    drift=[-1,*shares]
    basis=flavor["basis"]
    entries=((0,1),(1,2),(0,2),(1,1))
    terms=[]
    for signs in product(range(2),repeat=3):
        proj=flavor["identity"]
        for i in order:
            proj=proj*modes[i][signs[i]][1]
        vals=basis*proj*basis.transpose()
        lam=[modes[i][signs[i]][0] for i in range(3)]
        speed=sum(x*v for x,v in zip(drift,lam))
        terms.append((lam[0],speed,[vals[i,j] for i,j in entries]))
    radii=[sum(abs(speed)*abs(col[k]) for _,speed,col in terms) for k in range(4)]
    boxes=[(arb(1)/3 if start is None else start,end)]
    count=0
    while boxes:
        lo,hi=boxes.pop()
        mid=(lo+hi)/2
        width=(hi-lo)/2
        values=[flavor["radial"][i,j] for i,j in entries]
        for start,speed,col in terms:
            phase=acb(0,start+speed*mid).exp()
            values=[v+phase*c for v,c in zip(values,col)]
        data=[]
        for v,radius in zip(values,radii):
            err=abs(width)*radius
            data.append(v+acb(arb(0,err),arb(0,err)))
        us,cb,ub,cs=data
        a,b,u=[abs(v)**2 for v in data[:3]]
        q=us*cb*ub.conjugate()*cs.conjugate()
        readings=((a/(1-u)).sqrt(),(b/(1-u)).sqrt(),u.sqrt(),(q+a*b*u/(1-u)).arg())
        count+=1
        if any(abs(v-target)>d for v,target,d in zip(readings,flavor["shell"]["centers"],tolerances)):
            continue
        assert width>arb("1e-10"),(atoms,order,mid.str(12))
        boxes.extend(((lo,mid),(mid,hi)))
    return count


def verify():
    cost,strength,*_ = ordered["previous"]["setup"]()
    remaining = arb(1)/3
    cube = [c+arb(0,arb(2)**-28) for c in ordered["previous"]["box"]["centers"]]
    tolerances = [arb(".000055")]*3+[arb(".006")]
    data = {}
    def atom(n):
        if n not in data:
            height = acb.zeta_zero(n).imag
            x = arb(0) if n==center else study["center"]/(height**2+study["offset"])-1
            c = arb(0) if n==center else cost(x)[0]
            data[n] = x,c,study["native"]["plane_response"](height)
        return data[n]
    def pair_price(l,r):
        x,c,_ = atom(l)
        y,d,_ = atom(r)
        _,wl,wr = weights(x,y,remaining)
        return wl*c+wr*d
    winner = (11410,11432)
    least = pair_price(*winner)
    # 若最近的另一个编号也超过门槛，更远处由代价严格凸性排除。
    low,high = 11400,11500
    assert pair_price(low,center+1)>least and pair_price(center-1,high)>least
    successes=[]
    count=0
    checked_boxes=0
    for l in range(low+1,center):
        for r in range(center+1,high):
            price = pair_price(l,r)
            if (l,r)!=winner:
                assert 0 not in price-least
                if price>least:
                    break
            count+=1
            ns=(center,l,r)
            ws=weights(atom(l)[0],atom(r)[0],remaining)
            steps={n:ordered["evolution"](atom(n)[2],w,cube) for n,w in zip(ns,ws)}
            if (l,r)!=winner:
                max_mass=min(arb(1),(remaining*least/price).upper())
                modes=[spectral_modes(atom(n)[2],cube) for n in ns]
                shares=[w/remaining for w in ws[1:]]
                for order in permutations(range(3)):
                    checked_boxes+=exclusion(order,ns,modes,shares,max_mass,tolerances)
                continue
            for order in permutations(ns):
                obs=ordered["readings"](steps,order)
                errors=[abs(v-c) for v,c in zip(obs,flavor["shell"]["centers"])]
                if all(v<tol for v,tol in zip(errors,tolerances)):
                    successes.append((l,r,order))
                    print("Accepted",order,"weights",*[w.str(18) for w in ws],
                          "CKM over calibration box",*[o.str(12) for o in obs],flush=True)
                else:
                    assert any(v>tol for v,tol in zip(errors,tolerances)),(order,errors)
    assert successes==[(11410,11432,(11432,11431,11410))],successes
    print("Pairs at or below minimum",count,"continuous mass boxes",checked_boxes,"minimum port", (least*strength).str(21))
    assert least*strength<arb("1.5346448518982204e-29")


if __name__=="__main__":
    verify()
