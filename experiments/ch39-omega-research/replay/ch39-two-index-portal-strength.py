"""在相邻两枚 ζ 零点共同匹配引力后，中心电弱曲率需要多大的重门户。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat


def verify():
    root = Path(__file__).resolve().parent
    columns = run_path(str(root / "ch39-three-index-broken-gravity.py"))["gravity_columns"]
    n = 2_306_888_100_085
    max_cap = arb("3.427")
    mats,sites,heavy,heat = columns(root, (n,n+1,n-1), max_cap)
    r = (arb("246.22")/arb("1.22089e19"))**2/2
    g1,g2,yt = map(arb, (".467453",".63811",".8616"))
    light = ((-12,yt**2),(6,g2**2/2),(3,(g2**2+3*g1**2/5)/2))
    pi = arb.pi()
    for neighbour, bounds, unit_floor in (
        (1,("3.42428936","3.42428937"),"3.0326166e9"),
        (2,("3.42637240","3.42637241"),"9.3540244e9"),
    ):
        delta = mats[neighbour]-mats[0]
        mat = arb_mat([[1]*5] + [[a[i,j] for j in range(5)] for a in (mats[0],delta) for i in (0,1)])
        assert not mat.det().contains(0)
        weights = mat.solve(arb_mat([[1],[1],[arb("1e-122")],[0],[0]])).entries()
        def kernel(mass):
            return sum(w*heat(1,mass/x) for x,w in zip(sites,weights,strict=True))
        bare = -sum(d*k**2*kernel(k*r) for d,k in light)/(32*pi**2)
        coeff = [kernel(z+(-max_cap*r).union(max_cap*r))/(8*pi**2) for _,_,z in heavy]
        assert all(v > 0 for v in coeff)
        cost = sum(coeff)
        low, high = map(arb,bounds)
        assert bare-low**2*cost > arb(".258") > bare-high**2*cost
        assert bare-cost > arb(unit_floor)
        print("邻点相对编号",1 if neighbour==1 else -1,"单位门户下曲率的下界",
              (bare-cost).lower().str(16),"达到 .258 的最小最大门户区间",*bounds,
              "五权",*[w.str(12) for w in weights],flush=True)


if __name__ == "__main__":
    verify()
