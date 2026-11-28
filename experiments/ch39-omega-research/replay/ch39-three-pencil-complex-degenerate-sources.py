"""三项实壳铅笔重根内允许任意复相干时，六质量正向形状比的界。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb_mat, arb


def bounds(up, down, singleton, lower, upper):
    pair = [i for i in range(3) if i != singleton]
    for threshold, admissible in ((lower,False), (upper,True)):
        pencil = arb(threshold)*down-up
        block = acb_mat([[pencil[i,j] for j in pair] for i in pair])
        if admissible:
            assert pencil[singleton,singleton].real > 0
            assert block[0,0].real > 0 and block.det().real > 0
        else:
            assert block.det().real < 0


def verify(replay):
    orbit = run_path(str(replay/"ch39-two-zero-orbit-condition.py"))
    shell = run_path(str(replay/"ch39-shell-pencil-six-mass-window.py"))
    run_path(str(replay/"ch39-six-mass-extremal-shape.py"))
    model, up, down = orbit["calibration"]()
    frame, h = shell["shell_frame"](model)
    adj = model["adj"]
    a, b = h[0,1].real, h[1,2].real
    assert 0 < a < b
    k = (a*a+b*b).sqrt()
    root = arb(2).sqrt()
    # 列按 H 的本征值 -k、0、k 排序，首坐标均非零，故极点向量循环。
    spectral = acb_mat([[a/k/root,b/k,a/k/root], [-1/root,0,1/root],
                        [b/k/root,-a/k,b/k/root]])
    assert all(abs(spectral[0,j])>0 for j in range(3))
    sectors = [v*frame*spectral/abs(v.det())**(arb(1)/3) for v in (up,down)]
    gu, gd = [adj(v)*v for v in sectors]
    for j,(lo,hi) in enumerate((("34.28","34.31"),("3.43","3.45"),("1.16","1.17"))):
        bounds(gu,gd,j,lo,hi)
        print("H 的孤立线编号、允许任意复二阶块的最大迹比区间",j,lo,hi)
    # 纯极点 P 的二维空间是 Krylov 后两坐标，沿用全复子式。
    sectors = [v*frame/abs(v.det())**(arb(1)/3) for v in (up,down)]
    gu, gd = [adj(v)*v for v in sectors]
    positive = 34*gd-gu
    tail = acb_mat([[positive[i,j] for j in (1,2)] for i in (1,2)])
    assert positive[0,0].real>0 and tail[0,0].real>0 and tail.det().real>0


if __name__=="__main__":
    verify(Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent)
