"""给极点与三射线首站增加一站，区间核对五荷整个双钟域的轻带耦合。"""
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


ROOT = Path(__file__).resolve().parent


def certify(case):
    q, extra = case
    def lib(name):
        return run_path(str(ROOT/("ch39-"+name+".py")))
    graph = lib("three-ray-chiral-boundary")
    native = lib("native-light-projector-pump")
    lowest = lib("terminal-port-light-bands")["least_singular"]
    amplitude = lib("terminal-port-trial-energy")["restricted_lower"]
    cayley = lib("first-ray-rational-windows")["cayley"]
    adjoint = native["adjoint"]
    ctx.prec = 144
    sites = graph["geo"]["sites"]
    pole = graph["soft"]["pole"]
    ports = [pole, *graph["boundary_sites"]()]
    assert ports == [14,8,7,6]
    dist = [sum((x-y)**2 for x,y in zip(site,sites[pole],strict=True)) for site in sites]
    assert {i for i,v in enumerate(dist) if v == max(dist)} == set(range(6))
    rotations = [sites.index((y,z,x)) for x,y,z in sites]
    assert {rotations[i] for i in ports} == set(ports)
    assert {tuple(sorted({i,rotations[i],rotations[rotations[i]]})) for i in range(6)} == {(0,3,4),(1,2,5)}
    base, links = native["clock_matrices"](graph,q)
    sectors = []
    for frame in native["rotation_bases"](graph):
        restricted = acb_mat([frame.tolist()[i] for i in [*ports,extra]])
        sectors.append((restricted,adjoint(frame)*base*frame,
                        *(adjoint(frame)*m*frame for m in links)))

    @lru_cache(maxsize=8192)
    def cell(n,i):
        quadrant,k = divmod(i,n)
        mid = cayley(arb(2*k+1)/(2*n))
        rad = max(abs(mid-cayley(arb(t)/n)).upper() for t in (k,k+1))
        return acb(0,1)**quadrant * mid,rad

    pending = [(4,i,j) for i,j in product(range(16),repeat=2)]
    leaves = 0
    while pending:
        n,i,j = pending.pop()
        (u,ru),(v,rv) = cell(n,i),cell(n,j)
        move = arb(2).sqrt()*(ru+rv)+6*arb("1e-8")
        gaps, columns = [], []
        for restricted,fixed,a,b in sectors:
            w,gap = native["first_line"](fixed+a*u+adjoint(a*u)+b*v+adjoint(b*v))
            gaps.append(gap)
            columns.append((restricted*acb_mat([[z] for z in w])).entries())
        g = min(gaps)
        if g > 2*move:
            m = lowest(acb_mat(columns).transpose())
            if amplitude(m,move/(g-move)) > arb(".01"):
                leaves += 1
                continue
        assert n < 4096, (q,extra,n,i,j)
        pending.extend((2*n,2*i+a,2*j+b) for a,b in product(range(2),repeat=2))
    return q,extra,leaves


def actual_four_ports():
    def lib(name):
        return run_path(str(ROOT/("ch39-"+name+".py")))
    graph = lib("three-ray-chiral-boundary")
    forecast = lib("first-ray-forecast-phases")["FORECAST"]
    spectrum = lib("adjacent-pvm-transport")["spectrum"]
    ports = [graph["soft"]["pole"], *graph["boundary_sites"]()]
    for selected,n in forecast:
        for q,_ in forecast:
            levels,basis = spectrum(graph,q,n)
            assert levels[2] < 9 < levels[3]
            rows = acb_mat([[basis[i,j] for j in range(3)] for i in ports])
            values = (rows.conjugate().transpose()*rows).eig(algorithm="rump")
            assert all(x.is_finite() and x.imag.contains(0) and x.real > 0 for x in values)
            if q == selected:
                upper = min(x.real.sqrt().upper() for x in values)
                assert upper < arb(".000315")
                print("实际编号、指定荷、原四站最小奇异值上界",n,q,upper.str(12))
            else:
                assert min(x.real.sqrt().lower() for x in values) > arb(".03384")


if __name__ == "__main__":
    actual_four_ports()
    cases = list(product((1,-4,2,-3,6),(0,1)))
    with ProcessPoolExecutor(max_workers=len(cases)) as pool:
        for q,extra,count in pool.map(certify,cases):
            print("荷、附加站轨道代表、奇异值超过 0.01 的覆盖终格数",q,extra,count,flush=True)
