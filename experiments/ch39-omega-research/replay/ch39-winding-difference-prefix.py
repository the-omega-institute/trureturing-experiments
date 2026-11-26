"""连续实际零点的作用差、同高真空差与径向零位的辨识。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, arb_mat, ctx


root = Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
lib = run_path(str(root/"ch39-mincut-global-winding.py"))
shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
fixed = run_path(str(root/"ch39-height-independent-flavor-vacuum.py"))
angles = run_path(str(root/"ch39-mincut-reference-ckm-action.py"))
flavor = run_path(str(root/"ch39-chiral-neutral-flavor-phase.py"))


def rows():
    tri = lib["terms"](shell)
    target_x = arb("-.01046")
    target_y = target_x-angles["reference_angles"](flavor)[0]
    geometry = []
    equalities = []
    for q in range(1,106):
        equalities.append([arb(1), *fixed["gradient_columns"](lib,tri,q,target_x,target_y)])
        base = [arb(0),arb(0),arb(0)]
        for faces in tri.values():
            for z,j in faces:
                assert j in (0,1,2)
                v=1-(z**q).real
                base[0]+=v
                base[1]-=j*v
                base[2]+=int(j==2)*v
        unit=lib["coeff"](tri,q,arb(0))
        slope=[sum((z**q for z,j in tri[m] if j==1),acb(0)) for m in ((1,0),(0,1),(1,-1))]
        geometry.append((q,base,unit,slope))
    action,excitation,delta=[],[],[]
    for h in acb.zeta_zeros(1,103):
        h=h.imag
        phi=shell["geo"]["phases"](h)
        x,y=phi[724].arg(),phi[804].arg()
        s=arb(23)/25*arb(2)**-20/(h*h+shell["soft"]["point"]["study"]["offset"])
        full,gap=[],[]
        for q,(a,b,c),unit,slope in geometry:
            coeff=[v-s*z for v,z in zip(unit,slope,strict=True)]
            phase=lib["energy"](q,coeff,x,y)
            full.append(a+s*(b+s*c)+phase)
            gap.append(phase-lib["energy"](q,coeff,target_x,target_y))
        action.append(full)
        excitation.append(gap)
        delta.append(s)
    return list(map(list,zip(*equalities))),action,excitation,delta


def inverse_and_gain(equalities,reads,nuisance):
    k=len(nuisance[0])
    matrix=arb_mat([[*r,*([0]*k)] for r in equalities]+
                   [[*r,*s] for r,s in zip(reads,nuisance,strict=True)])
    inverse=matrix.inv()
    B=[[inverse[q,5+j] for j in range(len(reads))] for q in range(105)]
    gain=infinity_norm(B)
    return B,gain


def infinity_norm(rows):
    lengths=[sum(abs(a) for a in r) for r in rows]
    return max(v.lower() for v in lengths).union(max(v.upper() for v in lengths))


def difference_metrics(B):
    N=len(B[0])
    assert all(sum(r).contains(0) for r in B)
    # 参考广播的误差是独立测得的作用差的最大误差；共用能量误差另由完整行和给出。
    stars=[infinity_norm([[b for j,b in enumerate(r) if j!=anchor] for r in B]) for anchor in range(N)]
    neighbor=infinity_norm([[sum(r[k:]) for k in range(1,N)] for r in B])
    best=min(range(N),key=lambda k:float(stars[k]))
    assert all(stars[best] < v for k,v in enumerate(stars) if k!=best)
    return stars[0],best+1,stars[best],neighbor


def verify():
    ctx.prec=320
    equality,action,excitation,delta=rows()
    B,gain=inverse_and_gain(equality,action[:101],[[1]]*101)
    ref,best,bound,chain=difference_metrics(B)
    assert gain < arb("17.77837") and ref < arb("17.69430") and chain < arb("109.46186")
    assert best == 6 and bound < arb("17.24207") and gain/2 < arb("8.88919")
    print("原完整作用的共用零位、首枚差、最佳参考及相邻差",gain.str(14),ref.str(14),
          best,bound.str(14),chain.str(14),flush=True)
    _,direct=inverse_and_gain(equality,excitation[:100],[[] for _ in range(100)])
    B,gap_offset=inverse_and_gain(equality,excitation[:101],[[1]]*101)
    gap_ref,_,_,_=difference_metrics(B)
    assert direct < arb("12.77541") and gap_offset < arb("17.77837") and gap_ref < arb("17.69430")
    print("同高真空能差的直接读数、另有共用零位、相对首枚的能差",direct.str(14),gap_offset.str(14),gap_ref.str(14),flush=True)
    for k in (2,3):
        v=[(s/delta[0]) for s in delta[:100+k]]
        nuisance=[[t**j if j else arb(1) for j in range(k)] for t in v]
        _,amplitude=inverse_and_gain(equality,action[:100+k],nuisance)
        assert amplitude < arb("25.39000" if k==2 else "6.01518")
        print("独立径向多项式的系数个数、最短前缀、作用原始误差放大",k,100+k,amplitude.str(14),flush=True)


if __name__=="__main__":
    verify()
