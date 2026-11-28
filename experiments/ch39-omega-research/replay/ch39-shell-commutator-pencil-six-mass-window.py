"""壳循环、极点及有向交换子三项的任意实线性组合在六质量窗的谱源检验。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb


def verify(replay):
    pencil = run_path(str(Path(__file__).with_name("ch39-shell-pencil-six-mass-window.py")))
    orbit = run_path(str(replay/"ch39-two-zero-orbit-condition.py"))
    model, up, down = orbit["calibration"]()
    adj = model["adj"]
    frame, h = pencil["shell_frame"](model)
    gu, gd = [adj(z*frame)*(z*frame) for z in (up,down)]
    factor = abs(down.det()/up.det())**(arb(2)/3)
    k = 35*gd-factor*gu
    # 纯 P 时的重根允许整个二维复块。只有 P 和 i[H,P] 时另有孤立的第三线。
    block = acb_mat([[k[i,j] for j in (1,2)] for i in (1,2)])
    assert k[0,0].real > 0 and k[1,1].real > 0 and block.det().real > 0
    assert k[2,2].real > 0 and k[0,0].real*k[1,1].real > k[0,1].imag**2

    half = arb(1)/2
    cells = [(arb(-2),arb(2),0,False),(-half,half,0,True)]
    leaves = deepest = 0
    while cells:
        left, right, depth, reciprocal = cells.pop()
        t = left.union(right)
        f = pencil["eigenline"](t,h,reciprocal)
        # 在倒数坐标中另将虚方向除以 t，保持无穷端点非退化。
        if reciprocal:
            a, b, e = h[0,1].real, h[1,2].real, h[2,2].real
            w = (arb(0), a*(1-e*t), a*b*t)
        else:
            w = (arb(0), f[1], f[2])
        kf, kw = (pencil["form"](k,z) for z in (f,w))
        cross = f[0]*(k[0,1].imag*w[1]+k[0,2].imag*w[2])
        if kf>0 and kw>0 and kf*kw > cross**2:
            leaves += 1
            deepest = max(deepest,depth)
        else:
            assert depth < 40
            mid = (left+right)/2
            cells.extend(((left,mid,depth+1,reciprocal), (mid,right,depth+1,reciprocal)))

    f = pencil["eigenline"](arb(".9008065"), h, False)
    assert abs(f[0]) > 0
    multiplier = acb(1,arb(".007"))
    z = acb_mat([[f[0]], [f[1]*multiplier], [f[2]*multiplier]])
    assert factor*(adj(z)*gu*z)[0,0].real / (adj(z)*gd*z)[0,0].real > arb("34.3")
    shape = run_path(str(replay/"ch39-two-zero-source-correlations.py"))["shape"]
    masses = [m+arb(0,t.upper()) for m,t in zip(model["ref"][:6],model["tol"][:6],strict=True)]
    assert shape(masses[:3])/shape(masses[3:]) > arb("74.25")
    print("三项实铅笔的正谱源归一迹比上确界在 (34.3,35)，六质量要求大于 74.25。",
          "终格与最大深度", leaves, deepest)


if __name__=="__main__":
    verify(Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent)
