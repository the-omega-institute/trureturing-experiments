"""两个三高全程正解盒之间的输运：相对整体相位的最小改动秩恰为二。"""

from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    sixth_source = lib("adjacent-ordered-charge-six-cutoffs")["verify"](
        root,refine_node=2_306_888_100_086)[2][1]
    low, heavy, highest = [sixth_source[j].sqrt() for j in (0,-2,-1)]
    assert all(x is sixth_source[0] for x in sixth_source[:17]) and low<heavy<highest
    fifth = lib("fifth-third-pivot-three-high-positive")["endpoints"](root)[4]
    seventh = lib("fifth-third-pivot-em-colour-path")["solve"](
        root, lib("seventh-third-pivot-positive-path")["CENTERS"],(17,0,15,1,4,6),
        number=2_306_888_100_092)[4]
    highs = ((2,3,14),(1,4,6))
    shell = lib("three-ray-chiral-boundary")
    spectrum = lib("adjacent-pvm-transport")["spectrum"]
    first, last = [spectrum(shell,6,2_306_888_100_085+j)[1] for j in (5,7)]
    overlap = last.conjugate().transpose()*first
    ctx.prec = 192
    native_cos = sum(abs(overlap[i,j])**2 for i in highs[1] for j in highs[0]).sqrt()
    tails = []
    for boxes, bright in zip((fifth,seventh),highs,strict=True):
        assert [i for i,x in enumerate(boxes) if x>arb("1e-6")]==list(bright)
        leak = sum(x.sqrt()-low for i,x in enumerate(boxes) if i not in bright)
        tails.append((leak/(heavy-low)).upper().sqrt())
    # 重射线之间的角至少是原生高层间的角减去两个漏出角。
    # 逐步用余弦和平方根计算，检验每次相减后的角仍为正。
    cosine = native_cos.upper()
    for leak in tails:
        leak = leak.upper()
        assert cosine**2 + leak**2 < 1
        cosine = (cosine*(1-leak**2).sqrt() + (1-cosine**2).sqrt()*leak).upper()
    displacement = (2*(1-cosine)).sqrt()
    print("两组原生高层的重叠上界、两端漏出范数上界",
          native_cos.upper().str(14),*[v.upper().str(14) for v in tails],flush=True)
    print("重空间最大射线重叠上界、输运相对任意相位的第二奇异值下界",
          cosine.str(14),displacement.lower().str(14),flush=True)
    assert displacement > arb(".22")
    return displacement.lower()


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
