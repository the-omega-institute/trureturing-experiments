"""将两份双背景正解的共同分布分别冻结，在第三组真实零点检验。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

ROOT = Path(__file__).resolve().parent


def load(n):
    return run_path(str(ROOT/f"ch39-{n}.py"))


def verify():
    s = load("one-background-specific-line")
    old = load("shared-twelve-node-gauge-gap")
    shell = load("native-rh-ray-pvm")["shell_weights"]
    edges = old["paired_lines"](ROOT,shell,2_306_888_100_087)
    single,lib,(bg,) = load("two-zero-cutoff-gravity")["gravity_fields"](
        ROOT,(2_345_670,),tight_neutral=True)
    c = arb(7)/1000
    def point(i,t):
        a,b,w = edges[i]
        x = c*(a+t*(b-a))
        K,W = single["moments"](bg,lib["expint"],x)
        return arb_mat([[w*K],[w*W]])
    for labels,extra,ind,minima in s["layouts"]():
        shares = s["certify"](labels,extra,ind,minima)
        chosen = {i:point(i,t) for i,t in enumerate(labels[0]) if i!=ind}
        fixed = arb_mat([[-1],[-arb("1e-122")]])
        for v in chosen.values():
            fixed += v
        for coeff,(which,i,t) in zip(shares,extra,strict=True):
            if which is None:
                fixed += coeff*(point(i,t)-chosen[i])
        direction = (1,0) if ind==0 else (arb(1)/50,1)
        baseline = sum(x*y for x,y in zip(direction,fixed.entries(),strict=True))
        lows = []
        for j in range(256):
            t = (arb(j)/256).union(arb(j+1)/256)
            line = point(ind,t)
            lows.append(sum(x*y for x,y in zip(direction,line.entries(),strict=True)).lower())
        lower = baseline+min(lows)
        print("独立线、第三背景的最小引力残差投影下界",ind,lower.str(12),flush=True)
        assert lower > arb(".09" if ind==0 else ".0129")


if __name__ == "__main__":
    verify()
