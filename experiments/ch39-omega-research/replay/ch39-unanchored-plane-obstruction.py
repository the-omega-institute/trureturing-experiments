"""严格检验正电导两配置排除固定 Rayleigh 二平面，并验证变权标架。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb_mat, arb, ctx
import numpy as np


EDGES = ((1,6),(4,5),(5,7))


def verify(root):
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    light = run_path(str(root/"ch39-native-light-projector-pump.py"))
    frames = run_path(str(root/"ch39-third-edge-light-frame.py"))
    ctx.prec = 128
    a,(x,y) = light["clock_matrices"](shell,1)
    off = -x+y
    h = a+off+frames["adj"](off)
    columns = acb_mat(19,3)
    for k,(i,j) in enumerate(EDGES):
        columns[i,k] = 1
        columns[j,k] = h[j,i]
    update = columns * frames["adj"](columns)
    low, middle, large = arb(1)/100, arb(16), arb(256)
    cutoffs = []
    soft = 6*arb("1e-8")
    for weight in (low,middle,large):
        eigenvalues = [z.real for z in (h+(weight-1)*update).eig(multiple=True)]
        lo = sorted(z.lower() for z in eigenvalues)
        hi = sorted(z.upper() for z in eigenvalues)
        assert lo[3]-soft > hi[2]
        cutoffs.append((lo[1]-soft,hi[3],lo[2]-soft))
        print("等权、第二值下界、第三值上界、第四值下界及上界",weight.str(5),
              (lo[1]-soft).str(11),hi[2].str(11),(lo[3]-soft).str(11),hi[3].str(11))
    fraction = (middle-low)/(large-low)
    difference = cutoffs[1][0]-(1-fraction)*cutoffs[0][1]-fraction*cutoffs[2][1]
    assert 0 < fraction < 1
    assert difference > arb(".06962")
    assert cutoffs[2][2] > cutoffs[0][1]
    print("凸组合排除固定二平面的严格余量",difference.str(12))
    for charge in (1,-4,2,-3,6):
        base,links = light["clock_matrices"](shell,charge)
        nodes = np.exp(2j*np.pi*np.arange(8)/8)
        a,x,y = (frames["numeric"](v) for v in (base,*links))
        affine = np.stack([a + u*x+v*y + (u*x+v*y).conj().T for u in nodes for v in nodes])
        result = frames["verify_bands"](base,links,affine,EDGES)
        assert result is not None,charge
        print("荷、分段边、全部有理权区间的阈值整数和角格证书",charge,*result,flush=True)


if __name__=="__main__":
    root=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
    del sys.argv[1:]
    verify(root)
