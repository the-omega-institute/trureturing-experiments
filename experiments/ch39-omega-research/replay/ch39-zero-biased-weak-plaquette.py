"""用有向轮换读出右弱环，以正零点的原凯莱相位选取真空并反演高度。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx, fmpz_mat

flavour = run_path(str(Path(__file__).with_name("ch39-shell-plaquette-flavour.py")))
shell = flavour["shell"]


def capture(vertices, lengths, trees):
    lap, edges = shell["capture"](vertices,lengths,trees)
    return vertices, lap, edges


def geometry():
    _,_,_,(vertices,lap,edges) = shell["neutral"]["graph_constants"](port_builder=capture)
    neighbours,_ = flavour["cycles"](lap,edges)
    a,b = [i for i,p in enumerate(vertices) if p==(17,17,17)]
    mids = set(neighbours[a])
    assert mids == set(neighbours[b]) and len(mids)==3
    cycle = lambda p: p[1:]+p[:1]
    assert all(vertices.index(cycle(vertices[j])) in mids for j in mids)
    projected_edges = {tuple(sorted((vertices[i],vertices[j]))) for i,j in edges}
    assert len(set(vertices))==13 and len(projected_edges)==12
    # 连通图的顶点识别和重边合并保持连通。
    core = fmpz_mat(lap.nrows(),lap.ncols())
    for i in (a,b):
        for j in mids:
            core[i,i] += 1
            core[j,j] += 1
            core[i,j] -= 1
            core[j,i] -= 1
    def difference(i,j):
        return fmpz_mat([[int(v==i)-int(v==j)] for v in range(lap.nrows())])
    hub = difference(a,b)
    assert lap*hub==3*hub and core*hub==3*hub
    pairs = [(i,j) for i,j in combinations(range(lap.nrows()),2)
             if len(neighbours[i])==1 and neighbours[i]==neighbours[j]]
    assert len(pairs)==3
    for i,j in pairs:
        mode = difference(i,j)
        assert lap*mode==mode and core*mode==fmpz_mat(lap.nrows(),1)
    print("two hubs, three paths, three pairs of leaves",a,b,sorted(mids),pairs,flush=True)


def verify():
    geometry()
    ctx.prec = 160
    ratios = []
    for n in (11430,11431,11432):
        zero = acb.zeta_zero(n)
        height = zero.imag
        assert zero.real==arb(1)/2 and height>arb(1)/2
        bias = height/(height**2+arb(1)/4)
        assert ((zero-1)/zero).imag.overlaps(bias)
        r = 1/(1-arb(3).sqrt()*bias/2)
        ratios.append((3*r).sqrt())
        print(n,"height",height.str(13),"bias",bias.str(10),"mass ratio",ratios[-1].str(15),flush=True)
    assert acb.zeta_zero(1).imag>arb(1)/2
    a,b,c=ratios
    error = arb("1e-9")
    assert c < b*(1-error)/(1+error) and b*(1+error)/(1-error) < a
    print("neighbour gaps in mass ratio (ppb)", ((a/b-1)*10**9).str(10), ((1-c/b)*10**9).str(10),flush=True)


if __name__=="__main__":
    verify()
