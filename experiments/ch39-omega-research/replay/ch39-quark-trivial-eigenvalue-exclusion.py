"""用铅笔残差及下型正定性检验上型四站平凡二阶块的谱。"""
from pathlib import Path
from runpy import run_path

from flint import arb


def excluded(blocks, chart, lo, hi, masses, spectral_box):
    du,su,zu,ku,wu,_ = blocks[0]
    dd,sd,zd,kd,wd,_ = blocks[1]
    ratio = lo.union(hi)*masses[chart][-1]/masses[1-chart][-1]
    if chart == 0 and hi <= 0 and spectral_box[2] > masses[0][1] and su+ratio*sd > masses[0][1]:
        # 两个下型字符非负，使上型两个字符和高平凡本征值都高于中质量。
        return True
    sq = ku**2

    def constraints(low, high):
        r = (du-low)*(high-du)
        # C=1-uv 时，C R=|ku|²，且 R(u-v)wu=A。
        a = r*(zu-su) + sq*(su-low-high+du)
        d = r*(dd*(zd-sd)-kd**2) + sq*dd*sd
        char_d = (r-sq)*sd**2 - r*wd**2
        terms = [du-low, high-du, r-sq, (r-sq)*su**2+a*su-r*wu**2]
        if chart == 1:
            terms.extend((d + ratio*dd*a, char_d - ratio*sd*a))
        elif lo >= 0:
            terms.extend((ratio*d+dd*a, ratio*char_d-sd*a))
        elif hi <= 0:
            terms.extend((-ratio*d-dd*a, -ratio*char_d+sd*a))
        return terms

    # 对两个谱端点分别为一次式。四角同一条件均严格为负时，整个矩形被排除。
    g0,g1,h0,h1 = spectral_box
    corners = [constraints(g,h) for g in (g0,g1) for h in (h0,h1)]
    return any(all(y < 0 for y in row) for row in zip(*corners,strict=True))


def verify(root, spectral_box):
    pencil = run_path(str(root / "ch39-triple-geometry-pencil.py"))
    model = run_path(str(root / "ch39-global-independent-weak-frame-separation.py"))
    lib,geom,V,mass,A,B,_ = pencil["setup"](root)
    for chart,(base,other) in enumerate(((A,B),(B,A))):
        for branch in range(3):
            stack = [(arb(-1),arb(1),0)]
            count = depth = 0
            while stack:
                lo,hi,n = stack.pop()
                count += 1
                depth = max(n,depth)
                t = lo.union(hi)
                lines = lib["lines"](base,other,lo,hi,(branch,),geom)
                if branch in lines:
                    raw = pencil["physical_blocks"](lines[branch],chart,t,V,mass,geom,lib,())
                    if raw:
                        blocks = raw[0]
                        if n>4:
                            tight = model["centred"](lib,geom,V,mass,A,B,chart,branch,lo,hi)
                            if tight:
                                blocks = [tuple(c[j].intersection(x) if j in (0,1,2,4) else x
                                                for j,x in enumerate(b)) for c,b in zip(tight,blocks,strict=True)]
                        if excluded(blocks,chart,lo,hi,mass,spectral_box):
                            continue
                assert n < 25, (chart,branch,lo,hi)
                mid=(lo+hi)/2
                stack.extend(((lo,mid,n+1),(mid,hi,n+1)))
            print("图、本征支、成功覆盖的区间数、最大深度",chart,branch,count,depth,flush=True)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent, tuple(map(arb, ("0", ".00000615", ".010", ".018"))))
