"""共同树边权下双倍荷回路输运对两个夸克谱的绝对约束。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def verify():
    ctx.prec = 224
    root = Path(__file__).parent
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    common = run_path(str(root/"ch39-all-common-shell-ranks.py"))
    tree_model = run_path(str(root/"ch39-tree-weighted-higgs.py"))
    positive = run_path(str(root/"ch39-three-mode-higgs-hierarchy.py"))["positive"]
    ctx.prec = 224
    sites = shell["geo"]["sites"]
    rotation = [sites.index((y,z,x)) for x,y,z in sites]
    tree = tree_model["spanning_tree"](rotation)
    edges = {(i,j) for i,j,_,_ in shell["geo"]["edges"]}
    bases = list(common["character_bases"](shell))
    tree_lap = acb_mat(19,19)
    for i,j in tree:
        tree_lap[i,i] += 1
        tree_lap[j,j] += 1
        tree_lap[i,j] = tree_lap[j,i] = -1
    for basis in bases:
        level, _, rest = tree_model["lowest"](tree_lap,basis)
        assert level < arb(".10615") and all(v>arb(".31577") for v in rest)
    assert max(tree_lap[i,i].real.unique_fmpz() for i in range(19)) == 4
    assert max(sum(i in edge for edge in edges-tree) for i in range(19)) == 17
    first_height = acb.zeta_zero(1).imag
    offset = shell["soft"]["point"]["study"]["offset"]
    assert first_height>14 and offset>=0 and arb(23)/25*arb(2)**-20/14**2 < arb("1e-8")
    eps = arb(1)/1000
    gap = arb(".31577")-8*arb("1e-8")-arb(".10615")-34*eps
    assert gap>arb(".17561") and 68*eps/gap < arb(".388")

    height = acb.zeta_zero(11431).imag
    raw = shell["laplacian"](1,height)
    phase = {14: acb(1)}
    for _ in range(18):
        for i,j in tree:
            if (i in phase) != (j in phase):
                child,parent = (j,i) if i in phase else (i,j)
                phase[child] = -raw[child,parent]/abs(raw[child,parent])*phase[parent]
    assert len(phase) == 19
    projectors = {}
    identity = acb_mat([[int(i==j) for j in range(19)] for i in range(19)])
    for q in (1,4,-2):
        lap = shell["laplacian"](q,height)
        mat = acb_mat(19,19)
        for i,j in edges:
            weight = (1 if (i,j) in tree else eps)
            mat[i,j] = lap[i,j]*weight
            mat[j,i] = mat[i,j].conjugate()
            mat[i,i] += abs(mat[i,j])
            mat[j,j] += abs(mat[i,j])
        projector = acb_mat(19,19)
        for basis in bases:
            level,v,rest = tree_model["lowest"](mat,basis)
            assert level<arb(".2") and all(x>arb(".2") for x in rest)
            gauged = acb_mat([[phase[i]**(-q)*v[i,0]] for i in range(19)])
            projector += gauged*gauged.conjugate().transpose()
        projectors[q] = projector
    bounds = [arb(".0017900"),arb(".0022266")]
    for delta,bound in zip((projectors[1].conjugate()-projectors[1],
                            projectors[4].conjugate()-projectors[-2]),bounds):
        positive(identity*bound+delta)
        positive(identity*bound-delta)
    assert (arb(".8573")-arb(".01410"))/sum(bounds)>arb("209.92")
    print("全部正零点的共同低模数、任意剖面的归一质量差上界、11431 的质量差上界与拟合所需归一",
          3,"0.388","0.0040166",">209.92")


if __name__ == "__main__":
    verify()
