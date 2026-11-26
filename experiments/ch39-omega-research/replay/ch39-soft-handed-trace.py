"""在同一空间取向与同一实零点双钟下，三阶谱矩最早读出软辐边手性。"""
from collections import Counter
from pathlib import Path
from runpy import run_path
import sys

from flint import arb, acb, ctx


NUMBER = 2_306_888_100_085


def verify(root):
    exact = run_path(str(root/"ch39-mincut-real-triangle-tomography.py"))
    native = run_path(str(root/"ch39-unlabeled-native-connection.py"))
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 384
    sites = shell["geo"]["sites"]
    distances = [[sum((a-b)**2 for a,b in zip(x,y)) for x in sites] for y in sites]
    pos, mirrored = (native["reconstruct"](distances, i) for i in (0, 1))
    pole = pos[2]
    assert pos[3] == {0,3,4} and mirrored[3] == {1,2,5}
    replaced = {ij: (pos[-1][ij][0], k, sense, soften)
                for ij, (_,k,sense,soften) in mirrored[-1].items()}
    neg = (*mirrored[:-1], replaced)
    # 软边端点集互为反射，带权度及边平方权对第一、二阶迹的贡献相同。
    reflection = [sites.index((y,x,z)) for x,y,z in sites]
    for (i,j), (_,_,_,soft) in pos[-1].items():
        assert neg[-1][tuple(sorted((reflection[i],reflection[j])))][-1] == soft
    counts = Counter()
    all_faces = exact["enumerate_faces"]()
    for face, (k, z) in all_faces.items():
        mpos, mneg = (len(set(face) & G[3]) if pole in face else 0 for G in (pos, neg))
        assert 0 <= mpos <= 1 and 0 <= mneg <= 1
        if k == (0,0):
            # 无钟回路的反方向给出完全相同的实部。
            z = min(z, exact["conjugate"](z))
        counts[k,z] += mpos-mneg
    signed = {label: weight for label, weight in counts.items() if weight}
    assert signed == {((1,0),(510,1)): 3, ((1,0),(510,-1)): -3,
                      ((0,1),(1173,-2)): 3, ((0,1),(1173,2)): -3}
    h = acb.zeta_zero(NUMBER).imag
    phi = native["clocks"](pos,h)
    delta = arb(23)/25*arb(2)**-20/(h*h+shell["soft"]["point"]["study"]["offset"])
    ux, uy = [exact["unit"](z) for z in ((510,1),(1173,-2))]
    print("真实径向亏额",delta.str(17))
    for q, floor, ceiling in ((1,".27215586483699",".27215586483700"), (2,"1.07753","1.07755"),
                               (3,"2.38329","2.38331"), (4,"4.13603","4.13605"),
                               (6,"8.67851","8.67854")):
        ratio = -36*((ux**q).imag*(phi[1]**q).imag+(uy**q).imag*(phi[2]**q).imag)
        assert arb(floor) < ratio < arb(ceiling)
        plus, minus = (native["matrix"](G,q,phi,delta) for G in (pos, neg))
        for power in (1,2,3):
            difference = plus**power-minus**power
            moment = sum(difference[i,i] for i in range(19))
            assert moment.imag.contains(0) and moment.real.rad() < arb("1e-70")
            if power < 3:
                assert moment.real.contains(0)
            else:
                assert (moment.real/delta-ratio).contains(0)
                if q == 1:
                    assert arb("6.69878896654e-31") < moment.real < arb("6.69878896655e-31")
        print("荷、三阶迹差除以亏额、实际三阶迹差", q, ratio.str(17), (delta*ratio).str(17))


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
