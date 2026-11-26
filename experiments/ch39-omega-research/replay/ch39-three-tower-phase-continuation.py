"""在一个钟角开盒和径向亏额区间内同时延拓三座保轻重谱的正电导解。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, arb_mat, ctx


def rownorm(matrix):
    return max(sum(abs(matrix[i,j]) for j in range(matrix.ncols())).upper() for i in range(matrix.nrows()))


def sectors(sites):
    p=[sites.index((y,z,x)) for x,y,z in sites]
    cycles=[]
    rest=set(range(19))
    while rest:
        i=min(rest)
        cyc=[i,p[i],p[p[i]]] if p[i]!=i else [i]
        rest.difference_update(cyc)
        cycles.append(cyc)
    result=[]
    for k in range(3):
        kept=[cyc for cyc in cycles if k==0 or len(cyc)==3]
        basis=acb_mat(19,len(kept))
        for j,cyc in enumerate(kept):
            for l,i in enumerate(cyc):
                basis[i,j]=acb(0,2*arb.pi()*k*l/3).exp()/arb(len(cyc)).sqrt()
        result.append(basis)
    return result


def eigensystem(matrix,bases):
    rows=[]
    for basis in bases:
        block=basis.conjugate().transpose()*matrix*basis
        values,vec=block.eig(right=True,algorithm="rump")
        assert all(z.imag.contains(0) for z in values)
        full=basis*vec
        for i,z in enumerate(values):
            x=[full[j,i] for j in range(19)]
            unit=(sum(abs(a)**2 for a in x)).sqrt()
            x=acb_mat([[v/unit] for v in x])
            gap=min(abs(z.real-v.real).lower() for j,v in enumerate(values) if j!=i)
            rows.append((z.real,x,gap))
    rows.sort(key=lambda r:float(r[0]))
    assert len(rows)==19 and all(a[0]<b[0] for a,b in zip(rows,rows[1:]))
    return rows


def rayleigh(x,form):
    z=(x.conjugate().transpose()*form*x)[0,0]
    assert z.imag.contains(0)
    return z.real


def data(groups,geom,phase,weights,q):
    L=acb_mat(19,19)
    axes=[acb_mat(19,19) for _ in range(2)]
    forms=[]
    for j,orbit in enumerate(groups):
        H=acb_mat(19,19)
        for i,k in orbit:
            u=phase[i,k]
            H[i,i]+=1
            H[k,k]+=1
            H[i,k]=-u
            H[k,i]=-u.conjugate()
            if geom[i,k] in (724,804):
                axis=0 if geom[i,k]==724 else 1
                z=-acb(0,q)*u*weights[j]
                axes[axis][i,k]=z
                axes[axis][k,i]=z.conjugate()
        forms.append(H)
        L+=weights[j]*H
    return L,axes,forms


def targets(folded,c,perm):
    order=[*perm,*range(3,19)]
    signs=[-1,-1,1]+[1]*16
    factors=[1]*3+[32]*16
    branches=[-1]+[1]*18
    ys=[];slopes=[];center_slopes=[]
    for i,s,k,b in zip(order,signs,factors,branches):
        m=folded[i]
        r=(c*c+4*s*k*m).sqrt()
        assert r>0
        ys.append((c+b*r)/2)
        slopes.append(b*s*k/r)
        center_slopes.append((1+b*c/r)/2)
    return ys,slopes,center_slopes


def verify(root,source):
    prev=run_path(str(source))
    sh=run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    local=run_path(str(root/"ch39-mincut-reference-ckm-action.py"))
    flavor=run_path(str(root/"ch39-chiral-neutral-flavor-phase.py"))
    ctx.prec=288
    h=acb.zeta_zero(prev["NUMBER"]).imag
    clocks=sh["geo"]["phases"](h)
    x,y=[clocks[k].arg() for k in (724,804)]
    rho=arb("1e-9")
    dmax=arb("1e-12")
    dn=arb(23)/25*arb(2)**-20/(h*h+sh["soft"]["point"]["study"]["offset"])
    assert 0<dn<dmax
    assert sh["soft"]["point"]["study"]["offset"]>0 and arb(23)/25*arb(2)**-20/1000**2<dmax
    c=15-dn/3
    cw=c+arb(0,dmax/3)
    gap=local["reference_angles"](flavor)[0]
    ystar=arb("-.01046")-gap
    assert abs(y-ystar)+rho<arb("9e-6") and abs(x-y-gap)+2*rho<arb("2e-7")
    print("实际零点的钟角中心、实际亏额、统一亏额上端",*[z.str(19) for z in (x,y,dn,dmax)])
    groups=prev["orbits"](sh["geo"])
    geom={(i,j):d for i,j,d,z in sh["geo"]["edges"]}
    sites=sh["geo"]["sites"]
    rotate=[sites.index((y,z,x)) for x,y,z in sites]
    bases=sectors(sites)
    for q in (2,4,6):
        reference=sh["laplacian"](q,h)
        original=[abs(reference[i,j]) for (i,j),*_ in groups]
        phases={(i,j):-reference[i,j]/abs(reference[i,j]) for g in groups for i,j in g}
        num,variable,perm=prev["DATA"][q]
        weights=[arb(v)/10**15 for v in num]
        native,naxes,_=data(groups,geom,phases,original,q)
        fitted,faxes,forms=data(groups,geom,phases,weights,q)
        # 几何旋转保持相位；若反序则取共轭。钟边保持原定向，三条软边共轨。
        for i,j in geom:
            e=(i,j)
            a,b=rotate[i],rotate[j]
            other=tuple(sorted((a,b)))
            assert geom[e]==geom[other]
            if geom[e] in (724,804):
                assert a<b
            assert (phases[e]-(phases[other] if a<b else phases[other].conjugate())).contains(0)
        radial=[g for g,w in enumerate(original) if w<1]
        assert len(radial)==1 and len(groups[radial[0]])==3
        rhomat=-forms[radial[0]]
        assert all(sh["soft"]["pole"] in e and (phases[e]-1).contains(0) for e in groups[radial[0]])
        original_rows=eigensystem(native,bases)
        fitted_rows=eigensystem(fitted,bases)
        # 钟边每项的相位差至多 qρ，亏额的三辐边扰动范数至多 4Δδ。
        slopes=[[rownorm(A) for A in axes] for axes in (naxes,faxes)]
        delta_n=sum(slopes[0])*rho+4*dmax
        delta_f=sum(slopes[1])*rho
        lev=[a for a,u,g in original_rows]
        united=[a+arb(0,delta_n) for a in lev]
        masses=[a*abs(a-c) for a in lev]
        index=sorted(range(19),key=lambda i:float(masses[i]))
        massbox=[united[i]*abs(united[i]-cw) for i in index]
        assert all(a<b for a,b in zip(massbox,massbox[1:])) and massbox[2]<9<10<massbox[3]
        vals,_,_=targets([masses[i] for i in index],c,perm)
        _,yslopes,cslopes=targets(massbox,cw,perm)
        assignment=[index[k] for k in (*perm,*range(3,19))]
        J=arb_mat([[rayleigh(u,forms[g]) for g in variable] for z,u,gap0 in fitted_rows])
        P=J.inv()
        diff=arb_mat([[a-b] for (a,u,g),b in zip(fitted_rows,vals)])
        offset=rownorm(P*diff)
        variation=arb_mat(19,1)
        param_columns=[]
        for ax in range(2):
            ds=[]
            for (z,u,gd),j,py in zip(fitted_rows,assignment,yslopes):
                on,vec,gn=original_rows[j]
                v=united[j]
                assert v<cw or v>cw
                mprime=(2*v-cw)*(1 if v>cw else -1)
                # 谱投影的同特征标约化预解算子只需要本特征标内的间隔。
                err_n=slopes[0][ax]*(q*rho+2*delta_n/(gn-2*delta_n))
                err_f=slopes[1][ax]*(q*rho+2*delta_f/(gd-2*delta_f))
                assert gn>2*delta_n and gd>2*delta_f
                ds.append([rayleigh(u,faxes[ax])+arb(0,err_f)-
                           py*mprime*(rayleigh(vec,naxes[ax])+arb(0,err_n))])
            param_columns.append([row[0] for row in ds])
            direction=P*arb_mat(ds)
            for j in range(19): variation[j,0]+=abs(direction[j,0])*rho
        # 参考软边的导数是负的三辐边星图拉普拉斯，算子范数为 4。
        dd=[]
        for j,py,pc in zip(assignment,yslopes,cslopes):
            v=united[j]
            _,vec,gn=original_rows[j]
            sign=1 if v>cw else -1
            dl=rayleigh(vec,rhomat)+arb(0,8*delta_n/(gn-2*delta_n))
            dm=sign*((2*v-cw)*dl+v/3)
            dd.append(pc/3-py*dm)
        param_columns.append(dd)
        direction=P*arb_mat([[z] for z in dd])
        for j in range(19): variation[j,0]+=abs(direction[j,0])*dmax
        residual=offset+rownorm(variation)
        inc=[0]*19
        clock_change=[0]*19
        for j in variable:
            for a,b in groups[j]:
                inc[a]+=1;inc[b]+=1
                if geom[a,b] in (724,804):
                    clock_change[a]+=1;clock_change[b]+=1
        radius=2*residual
        change=2*max(inc)*radius+delta_f
        # D_wF 的变化逐行使用同特征标谱隙；每轨道 ||H||≤6，钟轨道相位扰动范数≤3qρ。
        Jerr=arb_mat(19,19)
        for j,(_,_,gd) in enumerate(fitted_rows):
            assert gd>2*change
            for k,g in enumerate(variable):
                Jerr[j,k]=12*change/(gd-2*change)
                if any(geom[e] in (724,804) for e in groups[g]):
                    Jerr[j,k]+=3*q*rho
        Pabs=arb_mat([[abs(P[i,j]) for j in range(19)] for i in range(19)])
        contraction=rownorm(Pabs*Jerr)
        print("荷、最小特征标内谱隙、预乘全盒残差、权盒半宽、统一压缩率",
              q,*[z.str(12) for z in (min(g for _,_,g in fitted_rows),residual,radius,contraction)],flush=True)
        assert contraction<arb(".5")
        contrast=arb(0)
        for k,(w,ref) in enumerate(zip(weights,original)):
            actual=w+arb(0,radius) if k in variable else w
            reference_box=arb(1).union(1-dmax) if k==radial[0] else arb(1)
            assert actual>0
            contrast=contrast.max((actual/reference_box).max(reference_box/actual))
        assert contrast<arb("8.001")
        # 对已求出的权盒也包住 D_p F。一个非零三行子式保证解对三背景坐标的导数满秩。
        for j,(_,_,gd) in enumerate(fitted_rows):
            for ax in range(2):
                widen=q*max(clock_change)*radius+2*slopes[1][ax]*(2*max(inc)*radius)/(gd-2*change)
                param_columns[ax][j]+=arb(0,widen)
        trio=(15,16,17)
        minor=arb_mat([[column[j] for column in param_columns] for j in trio]).det()
        floor=abs(minor).lower()
        assert floor>{2:233,4:2177,6:3808}[q]
        print("全盒电导反差上界、背景雅可比非零子式的三行（从1起）与绝对值下界",
              contrast.upper().str(12),(16,17,18),floor.str(12))
    print("钟角半宽与相对下密度",rho.str(10),(rho*rho/arb.pi()**2).str(15))


if __name__=="__main__":
    root=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
    source=Path(sys.argv[2]) if len(sys.argv)>2 else root/"ch39-three-even-positive-towers.py"
    verify(root,source)
