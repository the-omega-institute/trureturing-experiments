"""由固定正交补、行范数和共同酉的约束，排除第四壳的十六个混合三线轻上质量窗。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


def adj(matrix):
    return matrix.conjugate().transpose()


def data(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    x,y,_,readout,_,_,_ = lib("fourth-shell-two-gravity-rules")["shell_data"](root)
    shell = lib("three-ray-chiral-boundary")
    spec = lib("adjacent-pvm-transport")["spectrum"]
    n = lib("rank-two-generator-up-quark")["N"]
    _,labels,_ = lib("common-embedding-central-right-ckm")["verify_nodes"](
        root,((n+1,".0294"),),".1",(),((6,7,10,13),(6,7,10,12)))
    ctx.prec = 384
    mid = spec(shell,6,n+3)[1]
    last = adj(readout)
    sites = shell["geo"]["sites"]
    cycle = [sites.index((y,z,x)) for x,y,z in sites]
    groups = lib("charge-six-projector-algebra-transport")["groups"]
    trivial = groups(mid,cycle)[0]
    assert trivial==groups(last,cycle)[0] and len(trivial)==7
    basis = acb_mat([[(mid if j in trivial else last)[i,j] for j in range(19)] for i in range(19)])
    assert (adj(basis)*basis).contains(acb_mat([[int(i==j) for j in range(19)] for i in range(19)]))
    four = (shell["soft"]["pole"],*shell["boundary_sites"]())
    source,dst = [spec(shell,-4,n+j)[1] for j in (1,4)]
    beam = adj(y*y)*adj(x*x)*source
    weights = labels[n+1,-4]
    assert all(w>0 for w in weights)
    rho = acb_mat([[beam[i,j]*weights[j] for j in range(19)] for i in range(19)])*adj(beam)
    g = adj(basis)*rho*basis
    cross = adj(dst)*rho*basis
    rows = adj(dst)*basis
    original = adj(dst)*rho*dst
    fourrows = acb_mat([[dst[i,j] for j in range(19)] for i in four])
    return g,cross,rows,original,fourrows


def spectrum(matrix):
    values,columns = matrix.eig(right=True,algorithm="rump")
    order = sorted(range(len(values)),key=lambda j:values[j].real.mid())
    roots = [values[j].real for j in order]
    assert all(a<b for a,b in zip(roots,roots[1:]))
    sizes = [sum(abs(columns[i,j])**2 for i in range(matrix.nrows())).sqrt() for j in order]
    vec = acb_mat([[columns[i,j]/r for j,r in zip(order,sizes,strict=True)] for i in range(matrix.nrows())])
    return roots,vec


def bound(sector,support):
    g,cross,co,original,four = sector
    window = arb(".000008")
    ids = [i-1 for i in support]
    gram = acb_mat([[g[i,j] for j in ids] for i in ids])
    u = acb_mat([[co[i,j] for j in ids] for i in range(19)])
    linked = acb_mat([[cross[i,j] for j in ids] for i in range(19)])
    fixed = linked-u*gram
    ev,v = spectrum(gram)
    assert ev[0]>0
    z = fixed*v
    # 每个单行允许任意同范数方向。零乘子等价于行空间正交补的投影。
    shifts = [arb(0)]+[arb(10)**(arb(-10)+arb(16*j)/59) for j in range(60)]
    shifts.extend(-ev[0].mid()*arb(s) for s in (".98",".96",".85",".6",".2"))
    fixed_gram = original-linked*adj(u)-u*adj(fixed)
    base = [fixed_gram[i,i].real for i in range(19)]
    floor = []
    for i,b in enumerate(base):
        r2 = sum(abs(u[i,j])**2 for j in range(3))
        floors = [b-l*r2-sum(abs(z[i,j])**2/(ev[j]+l) for j in range(3)) for l in shifts]
        floor.append(max(t.lower() for t in floors))
    raw = acb_mat([[four[i,j]*floor[j] for j in range(19)] for i in range(4)])*adj(four)
    w = acb_mat([[1,0,0,0],[0,*([1/arb(3).sqrt()]*3)]])
    block = w*raw*adj(w)
    rest_trace = raw.trace().real-block.trace().real
    rest_square = (raw*raw).trace().real-(block*block).trace().real
    rest_gap = (2*rest_square-rest_trace**2).sqrt()
    assert (rest_trace-rest_gap)/2>window
    lv,svec = spectrum(block)
    if lv[0]>window:
        return lv[0].lower(),False
    assert lv[0]<window<lv[1]
    frac = (window-lv[0])/(lv[1]-lv[0])
    weighted_beam = adj(svec)*w*four
    eta = [max(arb(0),(abs(weighted_beam[0,j])*(1-frac).sqrt()
                       -abs(weighted_beam[1,j])*frac.sqrt()).lower())**2 for j in range(19)]
    weighted_u = acb_mat([[u[i,j]*eta[i] for j in range(3)] for i in range(19)])
    tau = adj(u)*weighted_u
    ts,t = spectrum(tau)
    linear = adj(t)*adj(weighted_u)*fixed*v
    constant = sum(b*h for b,h in zip(base,eta,strict=True))
    row_floor = sum(b*h for b,h in zip(floor,eta,strict=True))
    # 预选有理乘子。变换后矩阵仍酉，故每一行、每一列的模方和恰为一。
    pairs = []
    tick = arb("1e-7")
    for n in range(1,20):
        a = arb(n)/100000
        pairs.append(([a+tick,arb("8e-7"),arb("-7e-7")],[0,-a,-a]))
    pairs.append(([tick,arb(".0000101"),arb("-.0000019")],[0,arb(".000056"),0]))
    best = lv[0].lower()
    for x,y in pairs:
        coeff = [[ts[i]*ev[j]+x[i]+y[j] for j in range(3)] for i in range(3)]
        if all(c>0 for row in coeff for c in row):
            lower = constant-sum(x)-sum(y)-sum(abs(linear[i,j])**2/coeff[i][j]
                                                for i in range(3) for j in range(3))
            best = max(best,(lv[0]+lower-row_floor).lower())
    assert best>window,(support,best)
    return best,True


def verify(sector):
    selected = [(i,16,19) for i in (1,4,10,13)]
    selected += [tuple(sorted((16,18,i))) for i in (2,3,5,6,8,9,11,12,14,15,17,19)]
    found = [bound(sector,support) for support in selected]
    for support,(lower,refined) in zip(selected,found,strict=True):
        print("支撑、轻上窗反证下界、使用共同酉约束",support,lower.str(14),refined,flush=True)
    assert sum(refined for _,refined in found)==13
    return found


if __name__=="__main__":
    verify(data(Path(__file__).resolve().parent))
