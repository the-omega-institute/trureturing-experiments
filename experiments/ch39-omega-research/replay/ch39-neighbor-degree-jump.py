"""在两指定编号的对称近邻确定共同可加谱作用的第一次次数跃升。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def verify():
    root = Path(__file__).resolve().parent
    read = run_path(str(root/"ch39-light-spectrum-zero-selector.py"))
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    pump = run_path(str(root/"ch39-native-light-projector-pump.py"))
    lowest = run_path(str(root/"ch39-global-light-action-degree.py"))["lowest"]
    ctx.prec = 224
    blocks = list(pump["rotation_bases"](shell))
    adj = pump["adjoint"]
    displacements = {}
    for q,center in ((-4,2345666),(2,964890)):
        zeros = acb.zeta_zeros(center-83,167)
        assert len(zeros) == 167 and all(a.imag < b.imag for a,b in zip(zeros,zeros[1:]))
        moments = []
        for zero in zeros:
            matrix = shell["laplacian"](q,zero.imag)
            energy = [lowest(adj(b)*matrix*b) for b in blocks]
            moments.append(read["chebyshev_moments"](energy,4))
        for i,row in enumerate(moments):
            if i != 83:
                displacements[q,i-83] = [v-c for v,c in zip(row,moments[83],strict=True)]
    assert displacements[-4,-1][0] > 0 and displacements[2,1][0] < 0
    # e²+e = (81/8)T₂(x)+45T₁(x)+常数，其中 x=2e/9−1。
    quadratic = {k:-arb(81)/8*d[1]-45*d[0] for k,d in displacements.items()}
    assert all(v > arb(".06713626850") for k,v in quadratic.items() if abs(k[1])<=82)
    assert quadratic[-4,-83] < arb("-.03094371376")
    first_four = ((-4,-83),(-4,-59),(-4,-52),(2,-58))
    radius = read["interior_radius"]([displacements[k][:3] for k in first_four])
    assert radius > arb(".00011998")
    quartic = [-d[1]-d[2]-d[3]/4 for d in displacements.values()]
    assert len(quartic) == 332 and all(g > arb(".00355684951") for g in quartic)
    print("半径1至82的二次间隙、半径83的三次单形面距和四次间隙",
          min(v.lower() for k,v in quadratic.items() if abs(k[1])<=82).str(13),
          radius.str(13), min(v.lower() for v in quartic).str(13))


if __name__ == "__main__":
    verify()
