"""用原钟边与软辐边的轨道传递三个非恒定算术 Gram 源。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    graph = lib("three-ray-chiral-boundary")
    native = lib("native-light-projector-pump")
    adj = native["adjoint"]
    ctx.prec = 192
    pole = graph["soft"]["pole"]
    outer = sorted(graph["soft"]["star"]["orbit"])
    rays = graph["boundary_sites"]()
    assert (pole,outer,rays) == (14,[0,3,4],[8,7,6])
    identity = acb_mat([[int(i==j) for j in range(19)] for i in range(19)])
    away = acb_mat(identity)
    away[pole,pole] = 0
    # D 是软辐边亏额的负导数；压掉极点后恰为软辐边三个末端的坐标投影。
    derivative = acb_mat(19,19)
    for i in outer:
        derivative[i,i] += 1
        derivative[pole,pole] += 1
        derivative[pole,i] -= 1
        derivative[i,pole] -= 1
    mask = away*derivative*away
    assert all(mask[i,j] == int(i==j and i in outer) for i in range(19) for j in range(19))
    gram = lib("native-rh-ray-pvm")["gram_cholesky"](root)
    inputs = acb_mat([row[1:] for row in gram[1:]]).transpose()
    sigma = lib("terminal-port-light-bands")["least_singular"](inputs)
    assert sigma > arb(".0000899007") and arb(".31895")*sigma > arb(".0000286738")
    for q in (1,-4,2,-3,6):
        _, clocks = native["clock_matrices"](graph,q)
        for distance,clock in zip((724,804),clocks,strict=True):
            selected = mask*(clock+adj(clock))
            assert all(selected[i,pole].is_zero() for i in range(19))
            block = acb_mat([[selected[i,j] for j in rays] for i in outer])
            for source in rays:
                edges = [(i,j) for i,j,d,_ in graph["geo"]["edges"] if d==distance and i in outer and j==source]
                assert len(edges)==1 and abs(selected[edges[0]]).contains(1)
            for rows in (block.tolist(),block.transpose().tolist()):
                assert all(sum(not z.is_zero() for z in row)==1 for row in rows)
            assert all((adj(block)*block)[i,j].contains(int(i==j)) for i in range(3) for j in range(3))
            print("荷、钟边、三条射线传递到软辐边端点的轨道",q,distance,outer)
    print("三源边界限制奇异值下界与轻带传递的统一下界",sigma.str(14),".0000286738")


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
