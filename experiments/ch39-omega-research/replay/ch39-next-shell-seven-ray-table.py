"""用 Arb 重建 N+2 的上型振幅外包与 N+3 的原生射线有理外包。"""
import argparse
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


SCALE = 10**16
FILENAME = "ch39-next-shell-seven-ray-bounds.tsv"


def endpoints(value):
    """返回整数端点；除以 SCALE 后给出包含输入球的闭区间。"""
    scaled = value*SCALE
    return (int(scaled.lower().floor().unique_fmpz()), int(scaled.upper().ceil().unique_fmpz()))


def build(root):
    def load(name):
        return run_path(str(root/name))
    rank = load("ch39-rank-two-generator-up-quark.py")
    full = load("ch39-one-generator-gravity-full-ckm.py")
    center,radius = full["verify"](root)
    n = rank["N"]
    _,labels,_ = load("ch39-common-embedding-central-right-ckm.py")["verify_nodes"](
        root, ((n+1,".0294"),), ".1", (), ((6,7,10,13),(6,7,10,12)))
    shell = load("ch39-three-ray-chiral-boundary.py")
    spec = load("ch39-adjacent-pvm-transport.py")["spectrum"]
    ctx.prec = 384
    source,trained,target = [spec(shell,-4,j)[1] for j in (n+1,n+2,n+3)]
    adj,eye = rank["adj"],rank["eye"](19)
    x0 = load("ch39-one-generator-gravity-up-strange-bottom.py")["cayley"](
        root,rank["eye"],"ch39-full-ckm-gravity-cayley.tsv")
    h = acb_mat(19,19)
    directions = full["DIRECTIONS"]
    assert len(set(directions)) == len(center) == 18
    for step,(i,j,part) in zip(center,directions,strict=True):
        a = acb(0,step) if part=="i" else acb(step)
        h[i,j] += a
        if i!=j:
            h[j,i] -= a.conjugate()
    x = (2*(eye-h).inv()-eye)*x0
    radii = [abs((z+arb(0,radius.upper()))-z).upper() for z in center]
    dx = 2*sum(r*r*(1 if i==j else 2) for r,(i,j,_) in zip(radii,directions,strict=True)).sqrt()
    w = labels[n+1,-4]
    nominal = rank["read"](source,trained,w,adj(x)*adj(x))
    error = 2*dx*max(v.upper() for v in w).sqrt()
    rows = [("scale",SCALE), ("shells",n+2,n+3), ("sites",14,8,7,6)]
    for j,v in enumerate(nominal,1):
        interval = (max(arb(0),v.sqrt()-error)**2).union((v.sqrt()+error)**2)
        rows.append(("weight",j,*endpoints(interval)))

    assert (shell["soft"]["pole"],*shell["boundary_sites"]()) == (14,8,7,6)
    omega = acb(-1,arb(3).sqrt())/2
    for j in range(19):
        a = target[14,j]
        b = sum(target[i,j] for i in (8,7,6))/arb(3).sqrt()
        probs = [abs(target[8,j]+o*target[7,j]+o**2*target[6,j])**2/3 for o in (omega,omega.conjugate())]
        sectors = [abs(a)**2+abs(b)**2,*probs]
        active, = [k for k,p in enumerate(sectors) if p>0]
        assert all(p.contains(0) for k,p in enumerate(sectors) if k!=active)
        if active==0:
            z = a*b.conjugate()
            row = [n for p in (abs(a)**2,abs(b)**2,z.real,z.imag) for n in endpoints(p)]
            rows.append(("ray",j+1,*row))
        else:
            rows.append(("character",active,j+1,*endpoints(probs[active-1])))
    return "".join("\t".join(map(str,row))+"\n" for row in rows)


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write",action="store_true",help="写出重建的整数外包表")
    root = Path(__file__).resolve().parent
    write = parser.parse_args().write
    table = build(root)
    path = root/FILENAME
    if write:
        path.write_text(table,encoding="utf-8")
    elif path.read_text(encoding="utf-8")!=table:
        raise ValueError("Arb 复算结果与整数外包表不同")
    print("已写出整数外包表" if write else "Arb 复算逐项匹配整数外包表")
