"""固定同一质量单位，对十枚共用分界 9 的零点及一片正密度相位区做两圈复算。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

ROOT = Path(__file__).parent
fixed = run_path(str(ROOT / "ch39-chiral-clock-unit-readout.py"))
chiral = fixed["chiral"]


def verify():
    model,_,_,initial,_ = chiral["source"]["mean"]["inputs"]()
    ctx.prec = 192
    read = chiral["source"]["spectral_supplier"](model)
    hits = []
    for n in chiral["EXPECTED"][1:]:
        rows = read(acb.zeta_zero(n).imag)
        layers = {q:chiral["ordered"]([h for h,z in rows if z==q]) for q in (0,*model["charges"])}
        assert all(layers[q][2]<9<layers[q][3] for q in model["charges"])
        L,alpha = ("29.040",".1184") if n in (11661,11934) else ("29.166",".1212")
        print("实际零点",n)
        result = fixed["solve"](model,initial,layers,False,724,L,alpha)
        if arb(".117")<result<arb(".119"):
            hits.append(n)
        else:
            assert arb(".1211")<result<arb(".1214")
        if n==11661:
            anchor = layers
    assert hits==[11661,11934]

    # 每个钟角可偏移 1e-6 弧度。十二条钟边的六乘三扰动最多两项/行、四项/列。
    geom = chiral["source"]["geometry"]
    edges = [(i,j) for i,j,d,_ in geom["edges"] if d in (724,804)]
    assert len(edges)==12
    assert max(sum(a==i for a,b in edges) for i in range(19))<=2
    assert max(sum(b==j for a,b in edges) for j in range(19))<=4
    assert set(i for i,j in edges).isdisjoint(j for i,j in edges)
    # 这里只声称充分高零点的密度。两次三辐边修正的范数合计 <2e-12。
    shift = arb("2e-12")
    offset = model["shell"]["soft"]["point"]["study"]["offset"]
    assert offset>0 and 4*arb(23)/25*arb(2)**-20/(10000**2+offset)<shift/2
    assert acb.zeta_zero(11661).imag>10000
    delta = arb("1e-6")
    thick = {}
    for q,levels in anchor.items():
        move = arb(8).sqrt()*q*delta+shift
        # 谱在 [0,24]，|c−15| 微小，h 对 λ、c 的 Lipschitz 常数分别为 34、24。
        error = 34*move+24*2*shift
        thick[q] = [h+arb(0,error) for h in levels]
    assert all(thick[q][2]<9<thick[q][3] for q in model["charges"])
    print("11661 相位邻域",delta.str(6))
    result = fixed["solve"](model,initial,thick,False,724,"29.040",".1184")
    assert arb(".117")<result<arb(".119")
    density = (delta/arb.pi())**2
    print("合格零点的自然下密度界",density.str(12))


if __name__=="__main__":
    verify()
