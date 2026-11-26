"""在实零点21、22之间认证一只只包住指定 Weyl 点的正电导长盒。"""
from pathlib import Path
from runpy import run_path
import sys
from flint import acb, arb, ctx


def verify(root):
    sh=run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec=160
    center=arb("80.928425691033686017")
    weights=(arb("0.01615054064436668415"),arb("11.16447955978410845"))
    wr=arb("1e-15")
    g=sh["geo"]
    # 原钟块每类列互相正交，每列范数 sqrt(2)；三软边另计。
    freqs=[3*(g["theta"](d)-g["theta"](244)) for d in (724,804)]
    assert arb(2).sqrt()*sum(freqs)+arb("1e-7")<4
    ends=[acb.zeta_zero(k).imag for k in (21,22)]
    assert ends[0]>79 and sh["soft"]["point"]["study"]["offset"]>0
    stack=[tuple(ends)]
    checked=0
    while stack:
        a,b=stack.pop()
        if max(abs(a-center),abs(b-center))<arb("1e-13"):
            continue
        h=(a+b)/2
        H=sh["laplacian"](-3,h)
        for (i,j),w in zip(((0,2),(5,14)),weights):
            H[i,i]+=w-1
            H[j,j]+=w-1
            H[i,j]*=w
            H[j,i]=H[i,j].conjugate()
        vals = sorted((v.real for v in H.eig(multiple=True)),key=float)
        # 三四值在极小隙处接近浮点精度，次序用区间严格检查。
        assert all(x<y for x,y in zip(vals[:4],vals[1:5]))
        if vals[3]-vals[2] > 2*(4*(b-a)/2+4*wr):
            checked+=1
            continue
        assert b-a>arb("1e-17")
        stack.extend(((a,h),(h,b)))
    print("第21、22枚两端之间，严格排除额外交叉的高度区间数",checked)
    print("两条边的十进制中心及各自半宽",*(w.str(21) for w in weights),wr.str(5))
    print("与原局部唯一性盒衔接的中央高度半宽",arb("1e-13"))


if __name__=="__main__":
    root=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
    del sys.argv[1:]
    verify(root)
