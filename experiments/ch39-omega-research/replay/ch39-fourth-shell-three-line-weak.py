"""对七平凡源线与十二非平凡目的线的全部混合三线支撑判定 N+4 弱差。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb


def adj(m):
    return m.conjugate().transpose()


def norm(m):
    return sum(abs(z)**2 for z in m.entries()).sqrt()


def verify_data(root,mid,incoming,readout,energies,weights,*,trivial_counts=(1,2)):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    lo,h,k = energies
    target = adj(readout)
    shell = lib("three-ray-chiral-boundary")
    sites = shell["geo"]["sites"]
    rotation = [sites.index((y,z,x)) for x,y,z in sites]
    groups = lib("charge-six-projector-algebra-transport")["groups"]
    triv = groups(mid,rotation)[0]
    destination = groups(target,rotation)
    nontriv = sorted(destination[1]+destination[2])
    assert 15 in destination[0]
    high = incoming*acb_mat([[(h-lo).sqrt(),0],[0,(k-lo).sqrt()]])
    old = readout*high
    p0,p1 = [(arb(s)/arb("1.22089e19"))**2 for s in ("91.1876","1.176397e14")]
    weak = lib("full-neutral-free-cutoff")["weak_column"]
    weak_low = lib("native-rh-ray-pvm")["weak_series"](56)
    whole,large = lib("outer-shared-distance-kernel")["weak_bands"](weak_low,weak,p0,p1)
    ceiling = sum(c*(large if j==15 else whole) for j,c in enumerate(weights))
    assert ceiling<0 and lo*lo>p0
    limit = 16*p0.sqrt()
    current = [lo+norm(acb_mat([[old[j,0],old[j,1]]]))**2 for j in range(19)]
    assert min(z.lower() for z in current)>limit
    negative = sum(c*weak(z*z,p0,p1) for c,z in zip(weights,current,strict=True))
    assert negative<0
    eye = acb_mat([[int(i==j) for j in range(3)] for i in range(3)])
    cases = {}
    for ntriv in trivial_counts:
        off,yes = [],[]
        for src in combinations(triv,ntriv):
            for dst in combinations(nontriv,3-ntriv):
                name = (tuple(i+1 for i in src),tuple(i+1 for i in dst))
                l = acb_mat([[*[mid[i,j] for j in src],*[target[i,j] for j in dst]] for i in range(19)])
                a = adj(l)*high
                gram = adj(a)*a
                determinant = gram.det().real
                assert determinant>0
                coords = readout*l
                u = acb_mat([[coords[15,j] for j in range(3)]])
                b = acb_mat([[old[15,0],old[15,1]]])-u*a
                v0 = -b*gram.inv()*adj(a)
                r,needed = norm(u),norm(v0)
                if needed>r:
                    # 最小奇异值 >= sqrt(det(A*A)/tr(A*A))。
                    sig = (determinant/gram.trace().real).sqrt()
                    floor = lo+(sig*(needed-r))**2
                    assert floor>limit,(name,floor)
                    off.append((name,floor.lower()))
                    continue
                assert r>needed,(name,r,needed)
                # 三行二列的带号二阶子式组成单位左零行，且与最小范数解正交。
                kernel = acb_mat([[a[1,0]*a[2,1]-a[2,0]*a[1,1],
                                   a[2,0]*a[0,1]-a[0,0]*a[2,1],
                                   a[0,0]*a[1,1]-a[1,0]*a[0,1]]])
                v = v0 + kernel*((r*r-needed*needed).sqrt()/norm(kernel))
                delta = u-v
                denominator = r*r-(u*adj(v))[0,0]
                assert not denominator.contains(0)
                phase = eye-adj(delta)*delta/denominator
                assert (adj(phase)*phase).contains(eye)
                assert (u*phase-v).contains(acb_mat(1,3))
                corrected = old+coords*(phase-eye)*a
                assert corrected[15,0].contains(0) and corrected[15,1].contains(0)
                # 第 16 层精确等于训练底，其余逐层取本次同一酉矩阵的严格区间。
                amps = [lo+abs(corrected[i,0])**2+abs(corrected[i,1])**2 for i in range(19)]
                other = min(f.lower() for i,f in enumerate(amps) if i!=15)
                assert other>limit,(name,other)
                positive = weights[15]*weak_low(lo**2/p0)
                positive += sum(weights[i]*weak(f*f,p0,p1) for i,f in enumerate(amps) if i!=15)
                assert positive>0,(name,positive)
                yes.append((name,positive.lower(),other))
        cases[ntriv] = off,yes
        print("平凡线数、总数、排除数、弱差可达数",ntriv,len(off)+len(yes),len(off),len(yes),flush=True)
        print("可达支撑（源端平凡线；目的端非平凡线）",[name for name,_,_ in yes],flush=True)
        if off:
            print("排除支撑的第16层统一振幅下界",min(v for _,v in off).str(15),flush=True)
        if yes:
            print("可达支撑正端弱差与其余十八层振幅的统一下界",
                  min(v for _,v,_ in yes).str(13),min(v for _,_,v in yes).str(13),flush=True)
    print("恒等步弱差与第16层大截止时任意步弱差上界",negative.str(14),ceiling.upper().str(13),flush=True)
    return cases


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    _,_,incoming,readout,energies,weights,_ = lib("fourth-shell-two-gravity-rules")["shell_data"](root)
    shell = lib("three-ray-chiral-boundary")
    n = lib("rank-two-generator-up-quark")["N"]
    mid = lib("adjacent-pvm-transport")["spectrum"](shell,6,n+3)[1]
    return verify_data(root,mid,incoming,readout,energies,weights)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
