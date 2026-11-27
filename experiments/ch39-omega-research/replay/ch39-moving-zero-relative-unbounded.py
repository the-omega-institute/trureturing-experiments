"""同一实际高度既用于制备又用于读出时，两荷相对条件数的局部盲点证书。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    graph = lib("three-ray-chiral-boundary")
    pump = lib("native-light-projector-pump")
    ctx.prec = 192
    adj = pump["adjoint"]
    gram = lib("native-rh-ray-pvm")["gram_cholesky"](root)
    ports = [graph["soft"]["pole"],*graph["boundary_sites"]()]
    columns = []
    for row in gram[1:]:
        v = acb_mat([[row[ports.index(i)] if i in ports else 0] for i in range(19)])
        for w in columns:
            v -= w*(adj(w)*v)[0,0]
        norm = (adj(v)*v)[0,0].real.sqrt()
        assert norm > 0
        columns.append(v/norm)
    sources = acb_mat([c.entries() for c in columns]).transpose()
    assert all((adj(sources)*sources)[i,j].contains(int(i==j)) for i in range(3) for j in range(3))
    frames = list(pump["rotation_bases"](graph))
    def center(q,angles):
        fixed,links = pump["clock_matrices"](graph,q)
        source_rows = []
        for frame in frames:
            block = adj(frame)*fixed*frame
            for link,angle in zip(links,angles,strict=True):
                term = adj(frame)*link*frame*acb(0,angle).exp()
                block += term+adj(term)
            v,gap = pump["first_line"](block)
            assert gap > 1
            w = frame*acb_mat([[z] for z in v])
            source_rows.append((adj(w)*sources).entries())
        M = acb_mat(source_rows)
        return lib("first-zero-relative-quark-spectrum")["singular_values"](M)

    # 由已证荷 −4、软化为零的盲点盒取一个共同原相位的四次根分支。
    _,_,a,b = next(case for case in lib("first-ray-light-zeros")["CASES"] if case[0]==-4)
    angles = list(map(arb,(a,b)))
    for q in (-4,2):
        scale = arb(q)/(-4)
        midpoint = center(q,[scale*t for t in angles])
        shift = 2*arb(2).sqrt()*abs(scale)*arb(".0001")
        # 各字符最低线与格心投影之差至多 shift/(gap-shift)，三者直和取最大。
        error = shift/(1-shift)
        if q == 2:
            assert midpoint[0]-error > arb(".01")
        else:
            assert midpoint[-1]-error > arb(".1")
        assert center(q,[arb(0),arb(0)])[0] > arb(".001")
        print("荷、盲点盒心三源子空间最小/最大奇异值、投影扰动上界",
              q,midpoint[0].str(10),midpoint[-1].str(10),error.str(10))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
