"""中心三十次共同振幅冻结后，对两个相邻编号的十九截止与三代夸克作留出检验。"""
from pathlib import Path
from runpy import run_path
from flint import acb, acb_mat, arb, arb_mat, ctx


def basis(point,nodes,slopes,product):
    # 分开求每个 Lagrange 分子，允许留出谱值与中央位置重合。
    factors=[point-z for z in nodes]
    cardinal=[product(factors[:i]+factors[i+1:])/s for i,s in enumerate(slopes)]
    chi=product(factors)
    return cardinal,[chi*point**j for j in range(12)]


def orbital_rows(lap,positions,nodes,slopes,base):
    vals,vec=lap.eig(right=True,algorithm="rump")
    assert all(z.imag.contains(0) for z in vals)
    count=31
    blocks=[acb_mat(4,4) for _ in range(count)]
    for i,z in enumerate(vals):
        t=(z.real-18)/8
        lo,hi=basis(t,nodes,slopes,base["product"])
        column=acb_mat([[vec[k,i]] for k in positions])
        norm=sum(abs(vec[k,i])**2 for k in range(19))
        projector=(column*column.conjugate().transpose())*(1/norm)
        for j,s in enumerate(lo+hi):
            blocks[j]+=projector*s
    columns=[base["four_coordinates"](a) for a in blocks]
    return arb_mat(columns[:19]).transpose().tolist(),arb_mat(columns[19:]).transpose().tolist()


def inner(x,y):
    return sum(a.conjugate()*b for a,b in zip(x,y,strict=True))


def unit(x):
    norm=inner(x,x).real.sqrt()
    return [a/norm for a in x]


def read_mass(lap,pos,frame,nodes,slopes,base,evaluate):
    vals,vec=lap.eig(right=True,algorithm="rump")
    assert all(v.imag.contains(0) for v in vals)
    vecs=[]
    args=[]
    for i,z in enumerate(vals):
        n=sum(abs(vec[k,i])**2 for k in range(19)).sqrt()
        restricted=acb_mat([[vec[k,i]/n] for k in pos])
        vecs.append((frame.conjugate().transpose()*restricted).entries())
        args.append(basis((z.real-18)/8,nodes,slopes,base["product"]))
    lo,hi=zip(*args,strict=True)
    weights=evaluate(list(lo),list(hi)).entries()
    dominant=max(range(19),key=lambda i:float((abs(weights[i])*inner(vecs[i],vecs[i]).real).mid()))
    axis=unit(vecs[dominant])
    pivot=min(range(3),key=lambda i:float(abs(axis[i]).mid()))
    x=unit([int(i==pivot)-v*axis[pivot].conjugate() for i,v in enumerate(axis)])
    y=[(axis[(i+1)%3]*x[(i+2)%3]-axis[(i+2)%3]*x[(i+1)%3]).conjugate() for i in range(3)]
    axes=(x,y,axis)
    small=[]
    for i,v in enumerate(vecs):
        if i==dominant:
            small.append([acb(0),acb(0),acb(inner(v,v).real.sqrt())])
        else:
            small.append([inner(col,v) for col in axes])
    def entry(i,j):
        factors=[v[i]*v[j].conjugate() for v in small]
        return sum(f*w for f,w in zip(factors,weights,strict=True))
    a,b,c=entry(0,0).real,entry(0,1),entry(1,1).real
    u,v,h=entry(0,2),entry(1,2),entry(2,2).real
    # 重线单独写出，轻态求值保留仅为 1/h 阶的重块 Schur 修正。
    def schur(t):
        delta=h-t
        return a-t-abs(u)**2/delta,b-u*v.conjugate()/delta,c-t-abs(v)**2/delta
    def det(t):
        p,r,q=schur(t)
        return p*q-abs(r)**2
    cen=(a+c)/2
    rad=(((a-c)/2)**2+abs(b)**2).sqrt()
    light=[]
    for guess,sign in ((cen-rad,1),(cen+rad,-1)):
        middle=guess.mid()
        wid=arb("1e-12")
        while not (sign*det(middle-wid)>0>sign*det(middle+wid)):
            wid*=2
            assert wid<arb(".01"),(middle,wid,guess)
        light.append((middle-wid).union(middle+wid))
    assert light[0]<light[1] and abs(h)>arb("1e10")
    assert all(not (h-t).contains(0) for t in light)
    heavy=a+c+h-light[0]-light[1]
    assert not heavy.overlaps(light[0]) and not heavy.overlaps(light[1])
    allpairs=[]
    for t in light:
        p,r,q=schur(t)
        col=[-r,p] if abs(p).upper()>abs(q).upper() else [q,-r.conjugate()]
        col.append(-(u.conjugate()*col[0]+v.conjugate()*col[1])/(h-t))
        allpairs.append((abs(t),col))
    div=(heavy-a)*(heavy-c)-abs(b)**2
    allpairs.append((abs(heavy),[((heavy-c)*u+b*v)/div,(b.conjugate()*u+(heavy-a)*v)/div,1]))
    allpairs.sort(key=lambda pair:float(pair[0].mid()))
    masses=[p for p,_ in allpairs]
    assert 0<masses[0]<masses[1]<masses[2]
    lines=[]
    for _,col in allpairs:
        mapped=[sum(col[j]*axes[j][i] for j in range(3)) for i in range(3)]
        lines.append(unit(mapped))
    return masses,lines


def verify(root):
    common=run_path(str(root/"ch39-signed-three-source-quark-pair.py"))
    base=run_path(str(root/"ch39-hill-quark-common-polynomial.py"))
    ray=run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    native=run_path(str(root/"ch39-native-rh-ray-pvm.py"))
    ckm=run_path(str(root/"ch39-shell-clock-ckm-blocks.py"))
    leakage=run_path(str(root/"ch39-three-source-amplitude-leakage.py"))
    ctx.prec=416
    n=2_306_888_100_085
    pos=[ray["soft"]["pole"],*ray["boundary_sites"]()]
    assert pos==[14,8,7,6]
    height=acb.zeta_zero(n).imag
    lap=ray["laplacian"]
    def spectrum(h):
        zs=lap(6,h).eig(algorithm="rump")
        assert all(z.imag.contains(0) for z in zs)
        ordered=sorted(((z.real-18)/8 for z in zs),key=lambda z:z.mid())
        assert all(a<b for a,b in zip(ordered[:-1],ordered[1:],strict=True))
        return ordered
    nodes=spectrum(height)
    slopes=[base["product"](x-y for j,y in enumerate(nodes) if j!=i) for i,x in enumerate(nodes)]
    h=arb("5.66600e-35").union(arb("5.66604e-35"))
    k=arb(".635595244").union(arb(".635595246"))
    s=arb("3.291030650e-40").union(arb("3.291030654e-40"))
    u=[((8*t+18)/(8*nodes[-1]+18))**6000 for t in nodes[:-1]]+[arb(1)]
    cut=arb_mat([[(h+k*v/(v+s)).sqrt()] for v in u])
    frame=common["adapted_frame"](native["gram_cholesky"](root))
    change=base["compression"](frame,common)
    amplitudes,up,down,_=leakage["target_amplitudes"](ckm,common)
    rhs=[]
    L,A=[],[]
    for q,block in zip((-4,2),amplitudes,strict=True):
        low,high=orbital_rows(lap(q,height),pos,nodes,slopes,base)
        L+=low
        A+=high
        cols=arb_mat([[z] for z in common["coordinates"](block)])
        rhs+=change.solve(cols,algorithm="precond").tolist()
    A,L,b=arb_mat(A),arb_mat(L),arb_mat(rhs)
    def evaluate(low,high):
        dual=A.transpose().solve(arb_mat(high).transpose(),algorithm="precond").transpose()
        return dual*b+(arb_mat(low)-dual*L)*cut
    central_lead=evaluate([[0]*19],[[0]*11+[1]])[0,0]
    assert arb("-1.13089e14")<central_lead<arb("-1.13088e14")
    for label in (n-1,n+1):
        held=acb.zeta_zero(label).imag
        evaluations=[basis(x,nodes,slopes,base["product"]) for x in spectrum(held)]
        low,high=zip(*evaluations,strict=True)
        cutoffs=[abs(x)**2 for x in evaluate(list(low),list(high)).entries()]
        assert all(x>0 and (x>1 or x<1) for x in cutoffs)
        over=[i for i,x in enumerate(cutoffs,1) if x>1]
        expected=[1,2,3,4,6,7,8,9,10,19] if label<n else [1,2,3,4,7,8,9,10,18,19]
        assert over==expected
        lower,upper=("2.04740e29","2.04741e29") if label<n else ("9.24822e29","9.24825e29")
        assert arb(lower)<cutoffs[0]<arb(upper)
        print(label,"越过 M_P 的层号",over)
        print("十九层平方截止",*[x.str(8) for x in cutoffs],flush=True)
        data=[]
        for q,reference in zip((-4,2),(up,down),strict=True):
            mass,lines=read_mass(lap(q,held),pos,frame,nodes,slopes,base,evaluate)
            print("荷、三代 Yukawa 及二/三代比",q,*[v.str(10) for v in mass],(mass[1]/mass[2]).str(10))
            ratio=mass[1]/mass[2]
            observed=reference[1]/reference[2]*arb(".999")/arb("1.001")
            assert observed/ratio>(arb("4.5e9") if q==-4 else arb("1.2e12"))
            data.append(lines)
        def mix(i,j):
            return sum(a.conjugate()*b for a,b in zip(data[0][i],data[1][j],strict=True))
        s13=abs(mix(0,2)); c13=(1-s13*s13).sqrt()
        jarl=(mix(0,1)*mix(1,2)*mix(0,2).conjugate()*mix(1,1).conjugate()).imag
        reference=ckm["centers"][2]+ckm["sine_radius"]
        assert s13-reference>(arb(".02033") if label<n else arb(".04870"))
        assert jarl>0

        print("CKM 三正弦、有符号 J",*[v.str(10) for v in
              (abs(mix(0,1))/c13,abs(mix(1,2))/c13,s13,jarl)],flush=True)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
