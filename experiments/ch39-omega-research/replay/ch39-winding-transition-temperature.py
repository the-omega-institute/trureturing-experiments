"""三条振子链的首次计数转变、长链稳定温度与临界长度尾部。"""
from math import comb
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

free = run_path(str(Path(__file__).with_name("ch39-free-winding-cost.py")))
ORDER = 32


def bulk(temperature):
    # log(sinh(y)/y) 的交错级数，sin 的 2m 次矩为 C(2m,m)/4^m。
    value = -temperature.log()
    for m in range(1,ORDER+1):
        coefficient = (-1)**(m+1)*arb(2*m).zeta()/(m*(arb.pi()*temperature)**(2*m))
        value += coefficient*comb(2*m,m)/4**m
    next_order = ORDER+1
    remainder = arb(2*next_order).zeta()/(next_order*(arb.pi()*temperature)**(2*next_order))
    return value, remainder


def verify():
    ctx.prec = 170
    p, q = 23226, 92905
    angle = free["absolute"]["angle"](11431)
    alpha = 2*arb.pi()/3
    relative = arb("1e-18")*(arb(10000)/42)**4
    unit = relative*arb("1e-17")
    def static(count, band):
        d = abs(count*angle-band*alpha)
        return free["vacuum"](d) + 6*relative*(1-d.cos())/count
    saving = (static(p,1)-static(q,4))/(unit*(q-p))
    # 在两个指定温度核查有限振子和之差的符号。端点修正因同温而消去。
    def difference(t):
        center, remainder = bulk(t)
        radius = remainder*(q+p)/(q-p)
        return 6*t*(center+(arb(q)/p).log()/(2*(q-p))+arb(0,radius.upper()))-saving
    low, high = arb(".95215281352362"), arb(".95215281352363")
    assert difference(low)>0 and difference(high)<0
    # 在较高端点排除其他全部更长计数。任何相邻计数的余项也包括在内。
    center, err = bulk(high)
    min_rate = 6*high*(center-(2*p+1)*err)
    maximum_saving = static(p,1)-free["vacuum"](arb(0))
    bound = p+maximum_saving/(unit*min_rate)
    width = arb(".00003")
    assert free["vacuum"](width)-free["vacuum"](arb(0))>maximum_saving
    assert bound < 120000 and 120000*angle+width<6*alpha
    contenders = []
    for band in (1,2,4,5):
        left = int(((band*alpha-width)/angle).lower().floor().unique_fmpz())-1
        right = int(((band*alpha+width)/angle).upper().ceil().unique_fmpz())+1
        for length in range(max(p+1,left),right+1):
            distance = abs(length*angle-band*alpha)
            if length==q:
                assert distance<width
                contenders.append(length)
            else:
                assert distance>width
    assert contenders==[q]
    # 体积项符号决定长链稳定性。临界时剩下精确的 3 log p 热指数。
    critical_low, critical_high = arb("1.07339423805030"),arb("1.07339423805031")
    v,e = bulk(critical_low)
    assert v-e>0
    v,e = bulk(critical_high)
    assert v+e<0
    # 小长度以原始频率逐项核对有限三角矩公式。
    t = arb(".96")
    v,e = bulk(t)
    edge = ((1/t).sinh()*t).log()
    for length in (37,64):
        direct = sum((2*((arb.pi()*j/(2*length)).sin()/t).sinh()).log() for j in range(1,length))
        estimate = length*v+arb(length).log()/2+t.log()-edge/2
        assert abs(direct-estimate) < (length+1)*e
    print("首次转变温度",low.str(17),high.str(17),"所有可能较长竞争者截止",bound.str(12),
          "长链临界温度",critical_low.str(17),critical_high.str(17))


if __name__ == "__main__":
    verify()
