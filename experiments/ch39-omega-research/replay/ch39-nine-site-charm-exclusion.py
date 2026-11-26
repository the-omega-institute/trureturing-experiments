"""九个站点也无法同时产生参考顶底比与粲顶层级。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import arb,ctx


def verify():
    root = Path(__file__).parent
    flavor = run_path(str(root/"ch39-minimum-flavor-support.py"))
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 224
    up, down = flavor["flavor_maps"](shell,11431)
    adjoint, orbit = flavor["adjoint"], flavor["orbit"]
    reference = arb(".8573")/arb(".01410")
    k = adjoint(down)*down - 3*adjoint(up)*up/reference**2
    sites = shell["geo"]["sites"]
    rotation = [sites.index((y,z,x)) for x,y,z in sites]
    # 七阶及更小的主块正定，已由八站证书的 18.83 界保证。
    # 对可能含两个负方向的九阶块另外找一个正定八阶主块。
    bad_eight = set.union(*(orbit(g,rotation) for g in flavor["WITNESSES"]))
    count = 0
    bounds = []
    for group in combinations(range(19),9):
        cyclic = orbit(group,rotation)
        if group != min(cyclic):
            continue
        assert any(tuple(k for k in group if k != v) not in bad_eight for v in group)
        det = flavor["submatrix"](k,group).det().real
        if det < 0:
            bound = flavor["charm_bound"](up,down,group,reference)
            assert bound > arb(".208")
            bounds.append(bound)
            count += len(cyclic)
        else:
            assert det > 0
    assert len(bounds) == 209
    print("九站点中通过顶底必要条件的支撑数、循环等价类、粲顶比统一下界",
          count,len(bounds),min(bounds).str(14))


if __name__ == "__main__":
    verify()
