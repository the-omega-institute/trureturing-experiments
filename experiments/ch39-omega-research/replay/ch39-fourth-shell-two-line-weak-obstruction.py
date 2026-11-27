"""N+3 平凡荷六谱线与 N+4 非平凡线所张 84 个平面上任意 U(2) 的弱差上界。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb


def adj(a):
    return a.conjugate().transpose()


def norm(a):
    return sum(abs(v)**2 for v in a.entries()).sqrt()


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    _,_,incoming,readout,(low,h,k),weights,_ = lib("fourth-shell-two-gravity-rules")["shell_data"](root)
    shell = lib("three-ray-chiral-boundary")
    n = lib("rank-two-generator-up-quark")["N"]
    mid = lib("adjacent-pvm-transport")["spectrum"](shell,6,n+3)[1]
    target = adj(readout)
    sites = shell["geo"]["sites"]
    rotate = [sites.index((b,c,a)) for a,b,c in sites]
    groups = lib("charge-six-projector-algebra-transport")["groups"]
    src = groups(mid,rotate)[0]
    dest = groups(target,rotate)
    other = sorted(dest[1]+dest[2])
    assert 15 in dest[0]
    high_image = incoming*acb_mat([[(h-low).sqrt(),0],[0,(k-low).sqrt()]])
    original = readout*high_image
    p0,p1 = [(arb(v)/arb("1.22089e19"))**2 for v in ("91.1876","1.176397e14")]
    assert low*low>p0
    low_weak = lib("native-rh-ray-pvm")["weak_series"]()
    high_weak = lib("full-neutral-free-cutoff")["weak_column"]
    whole,large = lib("outer-shared-distance-kernel")["weak_bands"](low_weak,high_weak,p0,p1)
    ceiling = sum(c*(large if i==15 else whole) for i,c in enumerate(weights))
    print("只强制第16层超过16倍电弱尺度时的弱差上界",ceiling.upper().str(13),flush=True)
    assert ceiling<0
    rows,unresolved = [],[]
    for i in src:
        for j in other:
            pair = acb_mat([[mid[d,i],target[d,j]] for d in range(19)])
            a = adj(pair)*high_image
            coord = readout*pair
            # 未扰动的正交补行 b 与平面内的新行 vA 相加，其中 ||v|| 等于 ||u||。
            determinant = abs(a.det())
            adjugate = acb_mat([[a[1,1],-a[0,1]],[-a[1,0],a[0,0]]])
            frob = norm(a)
            assert frob>0
            bounds = []
            for row in range(19):
                if row in other and row!=j:
                    # 不同目的谱线正交，非平凡字符与源的平凡线亦正交。
                    fixed = low+sum(abs(original[row,l])**2 for l in range(2))
                    bounds.append(fixed.lower())
                    continue
                if row==j:
                    # 这条谱线属于所选平面，正交补贡献精确为零。
                    bounds.append((low+(determinant/frob)**2).lower())
                    continue
                u = acb_mat([[coord[row,0],coord[row,1]]])
                b = acb_mat([[original[row,0],original[row,1]]])-u*a
                radius = norm(u)
                # (b+vA) adj(A)=b adj(A)+det(A)v，对于奇异 A 亦成立。
                det_gap = abs(norm(b*adjugate)-radius*determinant)/frob
                direct = (norm(b)-radius*frob).lower()
                distance = max(arb(0),det_gap.lower(),direct)
                bounds.append((low+distance**2).lower())
            bands = [high_weak(v**2,p0,p1) if v>16*p0.sqrt() else whole for v in bounds]
            upper = sum(c*b for c,b in zip(weights,bands,strict=True)).upper()
            if upper>=0:
                unresolved.append((i+1,j+1,upper))
            rows.append((i+1,j+1,bounds[15],upper,sum(v>16*p0.sqrt() for v in bounds)))
    assert len(rows)==84
    worst = max(rows,key=lambda row:row[3])
    softest = min(rows,key=lambda row:row[2])
    print("平凡源层、非平凡目的层", [i+1 for i in src],[i+1 for i in other],flush=True)
    print("84 平面的第16层统一振幅下界及所在平面",softest[2].str(14),softest[:2],flush=True)
    print("各平面能保证高能的层数范围",min(row[4] for row in rows),max(row[4] for row in rows),flush=True)
    print("全部平面的弱差上界及所在平面",worst[3].str(14),worst[:2],flush=True)
    assert softest[2]>arb("1.28276e-8") and softest[2]>16*p0.sqrt()
    assert not unresolved,unresolved
    return rows


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
