"""用一只实际零点弦钟读出两类铰链负载及其形状残差。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def verify():
    ctx.prec = 192
    root = Path(__file__).parent
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 192
    target = 4-2*arb(3).sqrt()
    start = 10855
    heights = acb.zeta_zeros(start,1238)
    assert len(heights)==1238 and all(x.imag<y.imag for x,y in zip(heights,heights[1:]))
    probabilities = [(1+shell["geo"]["phases"](zero.imag)[724].real)/2 for zero in heights]
    distances = [abs(p-target) for p in probabilities]
    ordered = sorted(range(len(heights)),key=lambda j:float(distances[j]))
    first,second = ordered[:2]
    assert distances[first]>0 and distances[first]<distances[second]
    assert all(distances[second]<distances[j] for j in ordered[2:])
    assert (start+first,start+second)==(11366,10914)
    assert arb(".00008655165")<distances[first]<arb(".00008655166")
    print("窗口前两名与最小概率误差、11431 的概率",start+first,start+second,
          distances[first].str(12),probabilities[11431-start].str(12))
    for denominator,expected in ((100,19),(1000,1),(10000,1)):
        eps = arb(1)/denominator
        assert all(x<eps or x>eps for x in distances)
        density = ((2*(target-eps)-1).acos()-(2*(target+eps)-1).acos())/arb.pi()
        assert sum(x<eps for x in distances) == expected
        print("概率误差、窗口内命中数、全谱自然密度",eps.str(6),
              sum(x<eps for x in distances),density.str(12))


if __name__ == "__main__":
    verify()
