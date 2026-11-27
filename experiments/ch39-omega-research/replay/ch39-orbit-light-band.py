"""五荷三维轻带在三种三站轨道上的全环面稳定读取。"""
from concurrent.futures import ProcessPoolExecutor
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def certify(case):
    charge, half, root = case
    graph = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    pump = run_path(str(root/"ch39-native-light-projector-pump.py"))
    grid = run_path(str(root/"ch39-terminal-port-trial-energy.py"))
    ctx.prec = 144
    sites = graph["geo"]["sites"]
    boundary = sorted(graph["boundary_sites"]())
    assert boundary == [6,7,8]
    assert {sites.index((y,z,x)) for x,y,z in (sites[i] for i in boundary)} == set(boundary)
    base, links = pump["clock_matrices"](graph,charge)
    cells = grid["circle_cells"]()
    adj = pump["adjoint"]
    orbits=[]
    for start in range(len(sites)):
        x,y,z=sites[start]
        orbit=sorted({start,sites.index((y,z,x)),sites.index((z,x,y))})
        if len(orbit)==3 and orbit[0]==start:
            orbits.append(orbit)
    pole = sites[graph["soft"]["pole"]]
    distances = [sum((a-b)**2 for a,b in zip(site,pole,strict=True)) for site in sites]
    farthest = {tuple(orbit) for orbit in orbits if distances[orbit[0]] == max(distances)}
    assert max(distances) == 544 and farthest == {(0,3,4), (1,2,5)}
    sectors=[]
    for frame in pump["rotation_bases"](graph):
        orbit_columns=[next(k for k in range(frame.ncols()) if abs(frame[orbit[0],k])>arb(".5"))
                       for orbit in orbits]
        transformed=[adj(frame)*m*frame for m in links]
        clocks=[[t*z+adj(t*z) for z,_ in cells] for t in transformed]
        sectors.append((adj(frame)*base*frame,*clocks,orbit_columns))
    minimum=[[10**6]*3 for _ in orbits]
    for i in range(half*len(cells)//2,(half+1)*len(cells)//2):
        for j in range(len(cells)):
            change=arb(2).sqrt()*(cells[i][1]+cells[j][1])+6*arb("1e-8")
            for k,(fixed,first,second,indices) in enumerate(sectors):
                vector, gap = pump["first_line"](fixed+first[i]+second[j])
                assert gap > 2*change
                tilt = change/(gap-change)
                for g,index in enumerate(indices):
                    measured=grid["restricted_lower"](abs(vector[index]),tilt)
                    minimum[g][k]=min(minimum[g][k],grid["lower_millionths"](measured))
    return charge,half,minimum,orbits


def verify(root):
    cases=[(charge,half,root) for charge,half in product((1,-4,2,-3,6),range(2))]
    with ProcessPoolExecutor(max_workers=len(cases)) as pool:
        rows=list(pool.map(certify,cases))
    expected={(0,3,4):318950,(1,2,5):319740,(11,12,18):124950}
    assert set(expected).issubset(map(tuple,rows[0][3]))
    for g,orbit in enumerate(rows[0][3]):
        lower=min(min(r[2][g]) for r in rows)
        if tuple(orbit) in expected:
            assert lower>expected[tuple(orbit)]
            print("轨道与五荷三字符最小限制下界（百万分之一）",orbit,lower,flush=True)


if __name__=="__main__":
    import sys
    root=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
    verify(root)
