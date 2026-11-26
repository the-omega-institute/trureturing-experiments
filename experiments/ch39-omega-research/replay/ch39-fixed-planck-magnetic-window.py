"""固定普朗克单位或最大质量时，认证 1238 枚零点的每一个电弱一圈交点。"""
from collections import Counter
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


ROOT = Path(__file__).parent
mean = run_path(str(ROOT/"ch39-neutral-mean-local-mass.py"))
geometry = run_path(str(ROOT/"ch39-arc-relative-determinant.py"))


def spectral_supplier(model):
    # 仅十二条边带钟，单位电导谱求根后用 Weyl 范数界加入三辐边的微小减量。
    powers = {}
    for q in model["charges"]:
        mat = acb_mat(19,19)
        clock = []
        for i,j,d,u in geometry["edges"]:
            mat[i,i]+=1
            mat[j,j]+=1
            value = -u**q
            if d in (724,804):
                clock.append((i,j,d,value))
            else:
                mat[i,j]=value
                mat[j,i]=value.conjugate()
        powers[q]=(mat,clock)
    unit = model["shell"]["shell_spectra"]((0,),zero_index=11431)[0]
    def read(height, *, clocks=None):
        assert height>10000
        # 原有中性 11431 谱的微小包围和任一窗内谱的差小于 1e-12。
        neutral=[v+arb(0,"1e-12") for v in unit]
        center=sum(neutral)/18
        rows=[(v*abs(v-center),0) for v in neutral]
        if clocks is None:
            clocks=geometry["phases"](height)
        for q,(base,clock) in powers.items():
            mat=acb_mat(base)
            phase={d:z**q for d,z in clocks.items()}
            for i,j,d,value in clock:
                mat[i,j]=value*phase[d]
                mat[j,i]=mat[i,j].conjugate()
            values = mat.eig(multiple=True)
            assert len(values)==19
            for eig in values:
                assert eig.imag.contains(0)
                v=eig.real-arb(0).union(arb("1e-12"))
                rows.append((v*abs(v-center),q))
        assert len(rows)==113 and all(v>0 for v,q in rows)
        return rows
    return read


def electroweak_roots(rows,unit,model,initial,weights):
    mz,pi=model["mz"],arb.pi()
    lower,upper=((arb(x)/mz).log() for x in ("1e13","1e17"))
    thresholds=[]
    for v,q in rows:
        mass=v*unit
        assert mass>mz
        if mass<arb("1e17"):
            thresholds.append(((mass/mz).log(),q))
        else:
            assert mass>arb("1e17")
    thresholds.sort(key=lambda row:float(row[0]))
    assert all(l<r or q==u for (l,q),(r,u) in zip(thresholds,thresholds[1:]))
    slopes=model["beta"][:]
    rhs=[2*pi*v for v in initial]
    previous=lower
    found=[]
    for boundary,q in [*thresholds,(upper,None)]:
        if boundary>lower:
            a=previous.max(lower)
            diff=slopes[0]-slopes[1]
            if diff.contains(0):
                assert not (rhs[0]-rhs[1]-diff*a.union(boundary)).contains(0)
            else:
                scale=(rhs[0]-rhs[1])/diff
                if a<scale<boundary:
                    common=(rhs[0]-slopes[0]*scale)/(2*pi)
                    alpha=1/(common+(slopes[2]*scale-rhs[2])/(2*pi)+initial[2])
                    assert common>0 and alpha>0
                    found.append((scale,alpha))
                else:
                    assert scale<a or scale>boundary
        if q is not None:
            for i in range(3):
                slopes[i]+=weights[q][i]
                rhs[i]+=weights[q][i]*boundary
        previous=boundary
    return found


def verify():
    model,_,_,initial,weights=mean["inputs"]()
    ctx.prec=128
    read=spectral_supplier(model)
    planck=arb("2.435e18")
    # spectral_supplier 加载旧模块时可能设置更高精度。
    ctx.prec=128
    zeros=acb.zeta_zeros(10855,1238)
    assert len(zeros)==1238
    counts=[Counter(),Counter()]
    ranges={}
    for n,zero in enumerate(zeros,10855):
        rows=read(zero.imag)
        hmax=arb(0)
        for v,_ in rows:
            hmax=hmax.max(v)
        for j,unit in enumerate((planck,planck/hmax)):
            found=electroweak_roots(rows,unit,model,initial,weights)
            counts[j][len(found)]+=1
            for k,(L,alpha) in enumerate(found):
                key=(j,k)
                ranges[key]=alpha if key not in ranges else ranges[key].union(alpha)
                if n==11431:
                    print("11431 标定、交点序号、尺度及反推强耦合",j,k,(model["mz"]*L.exp()).str(8),alpha.str(10))
    assert counts==[{1:1238},{3:1238}]
    print("每种标定按交点次序的全窗反推强耦合")
    for key,a in sorted(ranges.items()):
        print(key,a.lower().str(12),a.upper().str(12))
        assert a<arb(".117") or a>arb(".119")


if __name__=="__main__":
    verify()
