"""在下一零点同时实现上型三个中心时，对保留相干的比例作联合下界和构造。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, arb


def two_levels(weights, entries):
    a,b,z = [sum(w*v for w,v in zip(weights,col,strict=True)) for col in entries]
    gap = ((a-b)**2+4*abs(z)**2).sqrt()
    return (a+b-gap)/2,(a+b+gap)/2


def spectral_weights(a, diagonal, error, r):
    eig = (r*a+(1-r)*diagonal).eig(algorithm="rump")
    roots = sorted((z.real.sqrt() for z in eig),key=lambda x:x.mid())
    assert all(z.imag.contains(0) for z in eig) and roots[0]>error
    assert all(x<y for x,y in zip(roots,roots[1:]))
    return [((w-error)**2).union((w+error)**2) for w in roots]


def exclude_mixed(weights,entries,chars,up,charm):
    w = [z.lower() for z in weights]
    # 一条平凡线与一条字符线。在每一阶，把剩余投影和的最小本征值取统一下界。
    flags = []
    for k in range(1,7):
        flag = min(two_levels([int(i in js) for i in range(7)],entries)[0].lower()
                   for js in combinations(range(7),7-k))
        flags.append(max(arb(0),flag))
    s = arb(".01")
    threshold = up+s*charm
    for js in combinations(range(13),7):
        wt = [w[j] for j in js]
        wc = [w[j] for j in range(13) if j not in js]
        trivial = wt[0]+sum((b-a)*f for a,b,f in zip(wt[:-1],wt[1:],flags,strict=True))
        for curve in chars:
            char = sum(a*b for a,b in zip(wc,curve,strict=True))
            assert trivial+s*char>threshold and char+s*trivial>threshold
    # 两条轻线分别属于两个字符，各自只有六条非零目的谱线。
    for js in combinations(range(12),6):
        first = [w[j] for j in js]
        other = [w[j] for j in range(12) if j not in js]
        x,y = [sum(a*b for a,b in zip(ws,c,strict=True)) for ws,c in
               zip((first,other),chars,strict=True)]
        assert x+s*y>threshold and y+s*x>threshold
    print("混合分配数、双字符分配数",1716,924,flush=True)


def construct(weights,entries,target,up,charm,top,t_text=".28994755"):
    # 在七个平凡目的层的第二、三条线上施加一次二阶旋转。
    mix = arb(t_text)
    perm = (2,6,5,3,4,0,1)
    wt = [weights[j] for j in perm]
    wt[1],wt[2] = (1-mix)*wt[1]+mix*wt[2],mix*wt[1]+(1-mix)*wt[2]
    lower = two_levels([v.lower() for v in wt],entries)
    upper = two_levels([v.upper() for v in wt],entries)
    e1,e2 = [v.union(w) for v,w in zip(lower,upper,strict=True)]
    assert e1<up<e2<charm
    root3 = arb(3).sqrt()
    o = acb(-1,root3)/2
    chars = []
    for phase,w in zip((o,o.conjugate()),(weights[13:],weights[7:13]),strict=True):
        strength = [abs(target[8,j]+phase*target[7,j]+phase**2*target[6,j])**2/3
                    for j in range(19)]
        curve = sorted((v for v in strength if v>0),key=lambda z:z.mid())
        assert len(curve)==6 and all(x<y for x,y in zip(curve,curve[1:]))
        chars.append(sum(a*b for a,b in zip(w,curve,strict=True)))
    c_high,c_middle = chars
    assert charm<c_middle<top<c_high
    print("与三个中心严格交错的四站谱",*[v.str(13) for v in (e1,e2,c_middle,c_high)],flush=True)


def verify(root):
    model = run_path(str(root/"ch39-next-shell-partial-coherence-up-threshold.py"))
    a,diagonal,target,error = model["source_data"](root)
    entries,char_bounds = model["target_data"](target)
    chars = [values[:6] for values in char_bounds]
    up,charm,top = map(arb,(".00000615",".00311",".8616"))
    lower = spectral_weights(a,diagonal,error,arb(".862503"))
    trace_cut = run_path(str(root/"ch39-trivial-block-coherence-trace-cut.py"))["clipped_minimum"]
    bound,old,new = trace_cut([v.lower() for v in lower[:7]],entries,up+charm)
    assert bound>up
    print("平凡块迹截面的最小值下界与两种顶点数",bound.str(14),old,new,flush=True)
    exclude_mixed(lower,entries,chars,up,charm)
    construct(spectral_weights(a,diagonal,error,arb(".862505")),entries,target,up,charm,top)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
