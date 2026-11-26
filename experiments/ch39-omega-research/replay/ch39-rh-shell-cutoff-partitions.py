"""固定极点与三射线 RH 谱权，将十九层任意合并至给定八个截止位置。"""
from bisect import bisect_left
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def subset_sums(values):
    result=[0]
    for v in values:
        result += [previous+v for previous in result]
    return sorted(result)


def distance_from_subsets(weights, target):
    # 只将精确整数送入组合搜索，统一扣去全部量化误差。
    unit=10**14
    center=int((unit*target.mid()).floor().unique_fmpz())
    approx=[int((unit*w.mid()).floor().unique_fmpz()) for w in weights]
    errors=sum(abs(w-arb(a)/unit).upper() for a,w in zip(approx,weights,strict=True))
    left= subset_sums(approx[:9])
    right=subset_sums(approx[9:])
    best=min(abs(a+right[j]-center) for a in left
             for pos in [bisect_left(right,center-a)]
             for j in (pos-1,pos) if 0<=j<len(right))
    return arb(best)/unit-errors-abs(target-arb(center)/unit).upper()


def verify(root):
    ctx.prec=192
    p0,p1=[(arb(e)/arb("1.22089e19"))**2 for e in ("91.1876","1.176397e14")]
    light=run_path(str(root/"ch39-full-neutral-free-cutoff.py"))
    low=arb("-.048342336655422").union(arb("-.048342336655420"))
    # 这里比较所有高位置与共同无截止矩；误差按每个实际探针的总变差另计。
    common=light["weak_column"](arb(1),p0,p1)
    regular_error=arb("7.468e-36")/arb("4e-15")
    first=-common/(low-common)
    half_width=arb("1e-8")
    integrator=run_path(str(root/"ch39-neutral-mz-cutoff.py"))["low_weak_column"]
    r=(arb("246.22")/arb("91.1876"))**2/2
    g1,g2,yt=map(arb,(".467453",".63811",".8616"))
    light_masses=[r*m for m in (g2**2/2,(g2**2+3*g1**2/5)/2,yt**2)]
    endpoint=[integrator(*light_masses,p0,p1,x_over_p0=(1+sign*half_width)**2) for sign in (-1,1)]
    low_error=abs(endpoint[0].union(endpoint[1])-low).upper()
    assert arb("2.865473795467")<first<arb("2.865473795470")
    shell=run_path(str(root/"ch39-native-rh-ray-pvm.py"))
    gaps=[]
    residuals=[]
    radii=[]
    gram=shell["gram_cholesky"](root)
    for q,levels,forward,reflected in shell["shell_weights"](root):
        bandgap=min((y-x).lower() for x,y in zip(levels[:-1],levels[1:],strict=True))
        for tag,probes in (("原序",forward),("镜序",reflected)):
            values=[]
            for i,weights in enumerate(probes,1):
                gap=distance_from_subsets([v.real for v in weights],first)
                values.append(gap.lower())
                rho=gram[i][i]
                margin=gap*abs(low-common)-max(regular_error,low_error)/rho
                residuals.append(margin.lower())
                # |Δw(S)| ≤ 9η/[ρ(g−η)]，g 是最小相邻带隙。
                relative=(margin*rho/abs(low-common)).lower()
                radii.append((bandgap*relative/(9+relative)).lower())
            gaps+=values
            print("荷、射线顺序、三枚枢轴在首个截止权上的严格距离下界",q,tag,
                  *[v.str(10) for v in values],flush=True)
    print("首节点权与全部取向的距离下界",first.str(16),min(gaps).str(12))
    assert min(gaps)>arb("1.2067e-5")
    print("低能带内单列最大偏移、弱差绝对下界、可承受的壳算子扰动半径",
          low_error.str(12),min(residuals).str(12),min(radii).str(12))
    assert min(residuals)>0 and min(radii)>0


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
