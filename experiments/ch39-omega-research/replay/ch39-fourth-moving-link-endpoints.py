"""分别移动联络的起点与终点，检验三种荷输运能否来自共同生成元。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, ctx


def adj(matrix):
    return matrix.conjugate().transpose()


def polar_bridge(source, target, old_groups, new_groups):
    overlaps = adj(target)*source
    phases = acb_mat(19,19)
    for old, new in zip(old_groups, new_groups, strict=True):
        for i,j in zip(old,new,strict=True):
            entry = overlaps[j,i]
            assert abs(entry)>0
            phases[j,i] = entry/abs(entry)
    return target*phases*adj(source)


def steps(root,y):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    n = lib("rank-two-generator-up-quark")["N"]
    shell = lib("three-ray-chiral-boundary")
    sites = shell["geo"]["sites"]
    rotation = [sites.index((b,c,a)) for a,b,c in sites]
    spectrum = lib("adjacent-pvm-transport")["spectrum"]
    groups = lib("charge-six-projector-algebra-transport")["groups"]
    result = []
    ctx.prec = 384
    for q,previous in zip((-4,2,6),(adj(y**2),y,y**3),strict=True):
        frames = [spectrum(shell,q,n+j)[1] for j in (2,3,4)]
        chars = [groups(frame,rotation) for frame in frames]
        incoming,outgoing = [polar_bridge(frames[i],frames[i+1],chars[i],chars[i+1]) for i in (0,1)]
        result.append(outgoing*previous*adj(incoming))
    up,down,gravity = result
    for actual,wanted in ((up,adj(down**2)),(gravity,down**3)):
        score = abs((adj(actual)*wanted).trace()/19)**2
        print("与荷二生成元的归一化迹模方",score.str(15),flush=True)
        assert score < 1
    return result


if __name__=="__main__":
    root=Path(__file__).resolve().parent
    run_path(str(root/"ch39-fourth-arithmetic-bridge-holdout.py"))["verify"](root,charge_steps=steps)
