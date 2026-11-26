"""全部带荷质量本征态与十八个重超荷模在匹配尺度以下的正谱积分。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

cut = run_path(str(Path(__file__).with_name("ch39-shell-vector-cut-budget.py")))
width = cut["width"]
census = run_path(str(Path(__file__).with_name("ch39-shell-charged-channel-census.py")))


def verify():
    ctx.prec = 256
    charged_matrix,neutral,r = census["mass_matrices"]()
    eigen,vectors = charged_matrix.eig(right=True)
    order = sorted(range(19),key=lambda j:float(eigen[j].real))
    charged = [eigen[j].real.sqrt()/r for j in order]
    assert all(lo < hi for lo,hi in zip(charged,charged[1:]))
    modes = []
    for j in order:
        norm = sum(abs(vectors[i,j])**2 for i in range(19)).sqrt()
        modes.append([vectors[i,j]/norm for i in range(19)])
    levels = [arb(0)]
    for v in sorted(width["neutral"]["shell"]["shell_spectra"]((0,))[0],key=float):
        if levels[-1].overlaps(v):
            levels[-1] = levels[-1].union(v)
        else:
            levels.append(v)
    total, count = arb(0),0
    for i,a in enumerate(charged):
        for j,other in enumerate(charged[:i+1]):
            profile = [v.conjugate()*w for v,w in zip(modes[i],modes[j])]
            for k,level in enumerate(levels[1:],start=1):
                b = level.sqrt()/r
                threshold = a+other+b
                if threshold > 1:
                    continue
                assert threshold < 1 and (other+b)**2-a**2 > arb(".00070")
                weight = cut["weight_bound"](neutral,profile,levels,k).upper()
                value = cut["cut_integral"](a,b,128,other).upper()
                multiplicity = 1 if i == j else 2
                total += multiplicity*weight*value
                count += multiplicity
    print("Open ordered charged pairs counted per neutral eigenvalue cluster",count)
    print("All charged pairs and all heavy hypercharge modes: upper per g_D squared",total.upper().str(12))
    assert count == 83 and total < arb(".19304")
    assert 2*arb(".19304") < arb(".4304")
    return total.upper()


if __name__ == "__main__":
    verify()
