"""由前22个真实正零点的完整磁拉普拉斯，以五荷共用权重制备算术三源。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx


ROOT = Path(__file__).resolve().parent
CHARGES = (1,-4,2,-3,6)
DENOMINATOR = 2**48


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    graph = lib("three-ray-chiral-boundary")
    pump = lib("native-light-projector-pump")
    ctx.prec = 256
    heights = [v.imag for v in acb.zeta_zeros(1,22)]
    theta = graph["geo"]["theta"]
    frequencies = [theta(d)-theta(244) for d in (724,804)]
    offset = graph["soft"]["point"]["study"]["offset"]
    data = []
    for h in heights:
        delta = arb(23)/25*arb(2)**-20/(h*h+offset)
        assert 0 < delta < arb("1e-8")
        row = [1,-delta]
        for q in CHARGES:
            for frequency in frequencies:
                phase = acb(0,q*h*frequency).exp()
                row.extend((phase.real,phase.imag))
        data.append(row)
    design = arb_mat(data)
    determinant = design.det()
    assert arb("1.69730") < determinant < arb("1.69731")
    inverse = design.inv()
    # 统一取软辐导数 D，或对各荷统一取第二组钟边的实部系数。
    selectors = ([1],[4+4*i for i in range(5)])
    rational = []
    for selected,limit in zip(selectors,("1.66701e9","2.862"),strict=True):
        weights = [sum(inverse[i,j] for i in selected) for j in range(22)]
        assert sum(abs(w) for w in weights) < arb(limit)
        rounded = [int((w*DENOMINATOR+arb(1)/2).floor().unique_fmpz()) for w in weights]
        fixed = [arb(n)/DENOMINATOR for n in rounded]
        assert all(abs(x-y) < arb(1)/(2*DENOMINATOR) for x,y in zip(weights,fixed,strict=True))
        rational.append(fixed)
        print("以 2^48 为分母的22个有符号整数",rounded)
    pole = graph["soft"]["pole"]
    away = acb_mat([[int(i==j and i!=pole) for j in range(19)] for i in range(19)])
    outer = set(graph["soft"]["star"]["orbit"])
    mask = acb_mat([[int(i==j and i in outer) for j in range(19)] for i in range(19)])
    gram = lib("native-rh-ray-pvm")["gram_cholesky"](root)
    ports = [pole,*graph["boundary_sites"]()]
    source = acb_mat([[gram[j][ports.index(i)] if i in ports else 0 for j in (1,2,3)]
                      for i in range(19)])
    assert sum(abs(z)**2 for z in source.entries()) < 4
    for q in CHARGES:
        laplacians = [graph["laplacian"](q,h) for h in heights]
        matrices = [sum((lap*w for lap,w in zip(laplacians,weights,strict=True)),acb_mat(19,19))
                    for weights in rational]
        target = mask*laplacians[0]*source
        rows = acb_mat([target.tolist()[i] for i in sorted(outer)])
        first_sigma = lib("terminal-port-light-bands")["least_singular"](rows)
        estimated = away*matrices[0]*away*laplacians[0]*source
        first_error = sum(abs(z)**2 for z in (estimated-target).entries()).sqrt()
        assert first_sigma > arb(".00050798") and first_error < arb("1e-11")
        assert arb(".31895")*first_sigma-first_error > arb(".00016202")
        print("荷、首零点送达软辐三端的三源奇异值下界",q,first_sigma.str(12))
        _,links = pump["clock_matrices"](graph,q)
        ideal = mask*(links[1]+pump["adjoint"](links[1]))*source
        rational_image = away*matrices[0]*away*matrices[1]*source
        error = sum(abs(z)**2 for z in (rational_image-ideal).entries()).sqrt()
        assert error < arb("1e-11")
        print("荷、实际22矩阵用固定有理权重后与理想三源矩阵的 Frobenius 误差上界",q,error.upper().str(12))
    assert arb(".31895")*arb(".0000899007")-arb("1e-11") > arb(".0000286738")
    print("设计矩阵行列式、全双钟域的三源最小奇异值下界",determinant.str(16),".0000286738")


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
