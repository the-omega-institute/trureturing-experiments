"""实际中央或右邻上，以显式正宽度的九峰三荷响应同时重建味和引力。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path
from flint import acb,arb,arb_mat,arb_series,ctx


CASES=(
    (2_306_888_100_085,(0,9,19,24,28,16,35),
     (".0000897345862",".8786001552",".00017001027",".1196993270",".2618020349",
      ".9207176041",".2956881667",".002005364316",".04017306517")),
    (2_306_888_100_086,(0,5,19,30,36,16,35),
     (".01481372753",".02359693038",".000185068488",".5697513892",".0000300649031",
      ".9241204889",".2960987637",".002005364316",".04017306517")),
)


def values_and_jac(fun,point):
    values=arb_mat([[w] for w in fun(point)])
    cols=[]
    for j in range(len(point)):
        jet=[arb_series([v,int(i==j)]) for i,v in enumerate(point)]
        cols.append([w[1] for w in fun(jet)])
    return values,arb_mat(cols).transpose()


def setup_weights(nodes,picks):
    anchors=[nodes[i] for i in picks]
    assert all(not (x-y).contains(0) for x,y in combinations(nodes,2))
    weights=[]
    narrow=arb(2)**-200
    for i,t in enumerate(nodes):
        if i in picks:
            weights.append([arb(0)]+[arb(int(i==j)) for j in picks])
        else:
            target=[narrow/(t-c)**2 for c in anchors]
            denom=1+sum(target)
            weights.append([1/denom]+[g/denom for g in target])
    return weights


def solve(root,case,field_rows,heat,moments):
    n,picks,start=case
    sparse=run_path(str(root/"ch39-sparse-rational-full-ckm.py"))
    native,model,A,B,mass,nodes,coords=sparse["setup"](root,n)
    ctx.cap=2
    reading=run_path(str(root/"ch39-native-rh-ray-pvm.py"))
    highroot=run_path(str(root/"ch39-ordered-cross-gravity-roots.py"))["high_pair"]
    high_weak=run_path(str(root/"ch39-full-neutral-free-cutoff.py"))["weak_column"]
    _,six,pivot,_=next(t for t in reading["shell_weights"](root,n) if t[0]==6)
    ctx.prec=256
    weights=[z.real for z in pivot[0]]
    abc=(sum(weights[:17]),weights[17],weights[18])
    mz,mp=arb("91.1876"),arb("1.22089e19")
    p0,p1=(mz/mp)**2,(arb("1.176397e14")/mp)**2
    low=arb(1).union(2)*p0
    high,slack=highroot(abc,low,field_rows,heat,moments,with_slack=True)
    high_amp=[x.sqrt() for x in high]
    base=low.sqrt()
    tar_basis=native["target_basis"](A,arb(".05")*B)*model["diag"]([1,1,acb(0,-1)])
    targets=[]
    for mat,ms in zip((A,B),mass,strict=True):
        t=model["adj"](tar_basis)*mat*tar_basis*ms[-1]
        assert t[0,1].contains(0) and t[1,2].imag.contains(0)
        targets.append((t[0,0].real,t[1,1].real,t[2,2].real,t[0,2].real,t[0,2].imag,t[1,2].real))
    levels=nodes[0]+nodes[1]+six
    spectrum=setup_weights(levels,(*picks,55,56))
    flow=[]
    for q,sec in enumerate(coords):
        flow.append([[sum(c[j]*row[k] for c,row in zip(sec,spectrum[19*q:19*q+19],strict=True))
                      for k in range(10)] for j in range(6)])
    def flavor(z,background,highs):
        strengths=[background,*z[:7],*highs]
        out=[[sparse["dot"](row,strengths) for row in rows] for rows in flow]
        return sparse["block_invariants"](out,z[-2],z[-1],targets)
    def at_center(z):
        return flavor(z,base.mid(),[v.mid() for v in high_amp])
    def uniform(z):
        return flavor(z,base,high_amp)
    point=list(map(arb,start))
    for _ in range(5):
        F,J=values_and_jac(at_center,point)
        dx=J.solve(F)
        point=[(x-y).mid() for x,y in zip(point,dx.entries(),strict=True)]
    radius=arb("1e-11")
    box=[x+arb(0,radius.upper()) for x in point]
    F,J=values_and_jac(at_center,point)
    pre=J.inv()
    F,_=values_and_jac(uniform,point)
    _,J=values_and_jac(uniform,box)
    defect=arb_mat([[int(i==j) for j in range(9)] for i in range(9)])-pre*J
    frac=max((abs(res)+sum(abs(v).upper() for v in row)*radius).upper()/radius
             for res,row in zip((pre*F).entries(),defect.tolist(),strict=True))
    assert frac<1
    assert all(base<x<arb(".948") for x in box[:7]+high_amp)
    assert box[-2]>0 and box[-1]>0 and box[-2]*box[-1]<1
    # 同一 SU(2) 旋转的充分条件：两对目标 Bloch 向量 Gram 相等，
    # 且两个目标轴不共线。正交标架可扩展成保持方向的 SO(3) 旋转。
    target_vectors=[]
    u,v=box[-2:]
    for p,s,T,re,im,w in targets:
        d=(T-u*v*s-(u-v)*w)/(1-u*v)
        b=(1-u*v).sqrt()
        target_vectors.append([2*re/b,-2*im/b,p-d])
    cross=sparse["cross"](*target_vectors)
    assert sum(v*v for v in cross)>0

    # 第十八、十九层的峰值精确；低层偏移由响应凸组合中全部峰的权重控制。
    leak=max(sum(row[1:]).upper() for row in spectrum[38:55])
    assert leak<arb("1e-48")
    slopes=[sum(abs(a) for g,a,z in field_rows)/(2*arb.pi()),
            sum(abs(g) for g,a,z in field_rows)/(4*arb.pi())]
    thermal_error=2*leak*max(slopes)*sum(abs(x) for x in weights[:17])
    assert thermal_error<slack
    weak=reading["weak_series"]()
    remaining=abc[1]*high_weak(high[0],p0,p1)+abc[2]*high_weak(high[1],p0,p1)
    end=[]
    for t in (arb(1),arb(2)):
        loc=(t*p0).sqrt()
        physical=loc.union(loc+leak)
        col=weak(physical**2/p0)
        end.append(sum(w*col for w in weights[:17])+remaining)
    assert end[0]>0>end[1]
    print(n,"九峰的味面占半宽、热核扰动及准许量、两端弱差",
          frac.str(7),thermal_error.str(6),slack.str(6),*[v.str(9) for v in end],flush=True)


def verify(root):
    gravity=run_path(str(root/"ch39-three-index-free-atoms.py"))
    moments,heat,fields=gravity["broken_fields"](root,[c[0] for c in CASES])
    for case,field in zip(CASES,fields,strict=True):
        solve(root,case,field,heat,moments)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
