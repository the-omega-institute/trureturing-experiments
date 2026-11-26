"""四维传播矩阵的极点差质量与全部阈值开启后的一圈预算。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx, fmpq

root = Path(__file__).parent


def verify():
    ctx.prec = 128
    shell = run_path(str(root/"ch39-zero-biased-weak-plaquette.py"))
    _,_,_,right = shell["shell"]["neutral"]["graph_constants"](port_builder=shell["capture"])
    vertices,lap,edges = right
    hubs = [i for i,v in enumerate(vertices) if v==(17,17,17)]
    assert len(vertices)==14 and len(edges)==15 and len(hubs)==2
    degrees = [int(lap[j,j]) for j in range(len(vertices))]
    assert all(degrees[i]==3 for i in hubs)
    length = 23226
    theta = run_path(str(root/"ch39-first-arithmetic-plaquette-winding.py"))["angle"](11431)
    delta = 2*arb.pi()/3-length*theta
    e = 2*(1-delta.cos())/length
    l = arb(4).union(arb(6))
    radial = 1/(1+e/(2*l))
    # 三条链在同一对公共极点之间形成平行规范质量边。
    profile_norm = length-1-arb(length**2-1)*(1-delta.cos())/(3*length)
    weight = 3*radial**2*profile_norm
    ratio = (3*radial+2*weight).sqrt()
    unextended = (3*radial).sqrt()
    allowed = 3*radial*((1+arb("1e-9"))**2-1)/(2*weight)
    assert 373 < ratio < 374 and allowed < arb("4.31e-14")
    print("等动能时极点差/叶差质量比",ratio.str(12),"原值",unextended.str(12),
          "1 ppb 容差给出的动能比上界",allowed.str(10))

    beta = [fmpq(d-22,3)+(length-1 if j in hubs else 0) for j,d in enumerate(degrees)]
    positive = sum(max(b,0) for b in beta)
    assert positive==fmpq(139312,3) and sum(beta)==fmpq(139072,3)
    weak_for_decade = (8*arb.pi()**2/(arb(positive)*arb(10).log())).sqrt()
    factor = (8*arb.pi()**2/(arb(".1")**2*arb(positive))).exp()
    assert weak_for_decade < arb(".0272") and factor < arb("1.186")
    print("极点站点各一圈系数",str(beta[hubs[0]]),"站点系数总和",str(sum(beta)),
          "一数量级所需对角耦合上界",weak_for_decade.str(10),
          "g_D=0.1 时尺度倍数上界",factor.str(10))


if __name__ == "__main__":
    verify()
