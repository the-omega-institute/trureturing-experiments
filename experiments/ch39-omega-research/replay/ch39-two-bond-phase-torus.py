"""在 12、32 两个匹配的全部对称与反对称相位上排除两壳联合引力四式。"""
from itertools import combinations, permutations
from pathlib import Path
from runpy import run_path
from subprocess import run
from tempfile import TemporaryDirectory

from flint import acb, arb


def prepare(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    fixed = lib("common-u2-eighth")
    unit = fixed["Q"]
    floor, ceil, midpoint = (fixed[key] for key in ("floor","ceil","mid"))
    lines = [[unit,fixed["PERIOD"]],[ceil(levels[2])]]
    for v in (levels[0], levels[1]-levels[0], levels[2]-levels[0]):
        lines.append([floor(v),ceil(v)])
    sites = lib("three-ray-chiral-boundary")["geo"]["sites"]
    groups = {12: [], 32: []}
    for i,j in combinations(range(19),2):
        length = sum((a-b)**2 for a,b in zip(sites[i],sites[j],strict=True))
        if length in groups:
            groups[length].append((i,j))
    assert len({i for e in (*groups[12],*groups[32]) for i in e}) == 12
    # 坐标置换对两种边稳定子的作用不同：只有 32 边的稳定子会交换两端。
    actions = [tuple(sites.index(tuple(p[k] for k in perm)) for p in sites)
               for perm in permutations(range(3))]
    for length, edges in groups.items():
        a,b = edges[0]
        orbit = {(min(t[a],t[b]),max(t[a],t[b])) for t in actions}
        assert orbit == set(edges)
        assert any(t[a]==b and t[b]==a for t in actions) == (length==32)
    geom = lib("fifth-seventh-shortest-site-gates")["geometry"]
    moment_errors = lib("fifth-seventh-cutoff-separation")["moment_errors"]
    for number in (2_306_888_100_092,2_306_888_100_093):
        base,_,_,basis,station,weights,heat,field = geom(root,fifth,columns,levels,number=number)
        target = moment_errors(field,heat,weights,[".75553"]*19)
        slopes = (sum(a for _,a,_ in field)/(2*arb.pi()),-sum(g for g,_,_ in field)/(8*arb.pi()))
        assert levels[2] < arb(".75553")
        for t, slope in zip(target, slopes, strict=True):
            tolerance = arb(".25")/abs(slope)
            lines.append([floor(t-tolerance), ceil(t+tolerance)])
        for w in weights:
            lines.append([floor(w),ceil(w)])
        for i in range(19):
            for j in range(2):
                coeff = []
                for edges in groups.values():
                    for sign in (1,-1):
                        coeff.append(sum((basis[a,i]+sign*basis[b,i]).conjugate()
                                         *(station[a,j]+sign*station[b,j])/2 for a,b in edges))
                terms = [base[i,j]-sum(coeff),*coeff]
                centers = [(midpoint(t.real),midpoint(t.imag)) for t in terms]
                rounding = sum(abs(t-acb(arb(a)/unit,arb(b)/unit))
                               for t,(a,b) in zip(terms,centers,strict=True))
                # 单位相位的两个坐标各至多一刻度误差；复投影截断再计两刻度。
                error = ceil(rounding + (2*sum(abs(t) for t in coeff)+2)/unit)
                lines.append([q for pair in centers for q in pair]+[ceil(abs(t)) for t in coeff]+[error])
    lines.extend(fixed["quarter_wave"]())
    return "\n".join(" ".join(map(str,line)) for line in lines)+"\n"


def experiment(root, fifth, columns, levels):
    data = prepare(root,fifth,columns,levels)
    with TemporaryDirectory() as temp:
        exe = Path(temp)/"cover"
        run(["c++","-std=c++17","-O3",str(root/"ch39-two-bond-phase-torus.cpp"),"-o",str(exe)],check=True)
        run([str(exe)],input=data,text=True,check=True)


def verify(root):
    source=run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    experiment(root,*source["verify"](root,return_operator="spectral"))


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
