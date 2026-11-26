"""在轻重谱之间平移重子三重性，核验选择律、连续规范提升及站点局域性。"""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


study = run_path(str(Path(__file__).with_name("ch39-three-leg-yukawa.py")))
# 左手 Weyl 多重态，依次给出色维数、弱维数、整数超荷 6Y 和轻场 X。
fields = {"Q":(3,2,1,0), "u":(3,1,-4,-1), "d":(3,1,2,1),
          "L":(1,2,-3,-1), "e":(1,1,6,2), "N":(1,1,0,0)}
pairs = (("Q","u",1),("Q","d",-1),("L","e",-1),("L","N",1))


def anomaly(color, weak, y, x):
    dim = color*weak
    return (int(color==3)*weak*x, int(weak==2)*color*x, dim*y*y*x,
            dim*y*x*x, dim*x**3, dim*x)


def verify_charges():
    for modulus in (2,3):
        offsets = [(a,b) for a,b in product(range(1,modulus), repeat=2) if (a+b)%modulus]
        assert offsets == ([] if modulus==2 else [(1,1),(2,2)])
    for left,right,higgs in pairs:
        x, y = fields[left][3],fields[right][3]
        assert (x+y+higgs)%3==0
        for a,b in product((0,1), repeat=2):
            assert ((x+y+higgs+a+b)%3==0)==(a==b==0)
    assert (2*(fields["L"][3]+1))%3==0 and (3*fields["Q"][3]+fields["L"][3])%3!=0
    # X=2Y−B+3L_tau，所有家族的离散剩余相同。
    rows=[]
    for family in range(3):
        for field,(color,weak,y,x) in fields.items():
            shift = ({"L":3, "e":-3, "N":-3}.get(field,0) if family==2 else 0)
            rows.append(anomaly(color,weak,y,x+shift))
    assert [sum(terms) for terms in zip(*rows)] == [0]*6
    for color,weak,y,x in fields.values():
        pair = (anomaly(color,weak,y,x+1), anomaly(color,weak,-y,-x-1))
        assert all(a+b==0 for a,b in zip(*pair))
    print("三代加单态中微子的六项连续反常为零，各重 Dirac 对逐对抵消")


def actual_projectors(number):
    h = acb.zeta_zero(number).imag
    lap = study["shell"]["laplacian"]
    neutral = lap(0,h)
    center = sum(neutral[i,i].real for i in range(19))/18
    floors = []
    for q in (1,-4,2,-3,6):
        v = study["edge"]["mode_columns"](study["edge"]["filters"]["spectra"](q,h,center))
        p = v*study["adjoint"](v)
        # 全部站点对检验，故与 P 对易的站点对角算子必须是标量。
        floors.append(min(abs(p[i,j]).lower() for i in range(19) for j in range(i)))
    assert min(floors)>arb("4.3e-6")
    print(number,"五种荷的轻投影最小非对角元绝对值下界",*(f.str(8) for f in floors))


def polynomial_degree():
    low, high = [], []
    for n in (22027,28325,11661,11934,12207):
        h = acb.zeta_zero(n).imag
        lap = study["shell"]["laplacian"]
        neutral = lap(0,h)
        center = sum(neutral[i,i].real for i in range(19))/18
        for q in (1,-4,2,-3,6):
            levels = [abs(r.real*(r.real-center))**2 for r in lap(q,h).eig(algorithm="rump")]
            assert sum(x<81 for x in levels)==3 and sum(x>81 for x in levels)==16
            low.extend(x for x in levels if x<81)
            high.extend(x for x in levels if x>81)
        if n in (28325,12207):
            for band in (low,high):
                vals = sorted(band,key=float)
                assert all(a<b for a,b in zip(vals,vals[1:]))
            print("共同原壳谱多项式：轻、重谱点及所需精确次数",len(low),len(high),len(low)+len(high)-1)


if __name__ == "__main__":
    ctx.prec = 448
    verify_charges()
    for number in (22027,28325):
        actual_projectors(number)
    polynomial_degree()
