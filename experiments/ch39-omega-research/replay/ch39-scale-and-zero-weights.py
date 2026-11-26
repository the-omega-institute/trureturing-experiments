"""四十一项算术读数对真实零点的共同尺度、逐点谱权与低处总谱权的限制。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, arb_poly, ctx


local = run_path(str(Path(__file__).with_name("ch39-three-mode-time-localization.py")))
clock = local["clock"]


def identify_pair(heights, windows, cutoff):
    (h1, window1), (h2, window2) = windows
    assert heights[8192]/h1*window1.lower().sqrt() > cutoff
    factor = h2/h1*(window1/window2).sqrt()
    assert all(lo < hi for lo,hi in zip(heights,heights[1:]))
    pairs = []
    j = 0
    for i,height in enumerate(heights[:8192]):
        target = height*factor
        assert target < heights[-1]
        while heights[j] < target:
            j += 1
        if j:
            assert heights[j-1] < target
        if heights[j].overlaps(target):
            pairs.append((i+1,j+1))
            assert target < heights[j+1]
        else:
            assert target < heights[j]
    assert pairs == [(1,2)]
    assert abs(window1.sqrt()-1) < arb("1.01e-16")
    print("尺度上限",cutoff.str(8),"允许的窗口编号对",pairs)
    return window1


def small_scale_weights(times, scale_squared, error):
    moments = clock["cleared_moments"](times,41)
    poles = local["partial_fractions"](41)
    q = arb_poly([1])
    for rate,_ in poles:
        q *= arb_poly([1,arb(rate)/200])

    def expectation(numerator):
        residues = local["residues_of"](numerator,poles)
        nominal = sum(c*t for c,t in zip(residues,times[1:],strict=True))
        return nominal+arb(0,error*sum(abs(c) for c in residues))

    nodes = [200/z.imag**2 for z in acb.zeta_zeros(1,13)]
    intervals = []
    for i,x in enumerate(nodes[:12]):
        p = local["window_polynomial"](moments,x)
        gaps = [(x-scale_squared*nodes[i+1]).lower()]
        if i:
            gaps.append((scale_squared*nodes[i-1]-x).lower())
        radius = min(gaps)*arb(".99")
        assert radius > 0
        squared = p**2
        minorant = squared*(1-(1/radius**2)*arb_poly([-x,1])**2)
        # 先平移多项式再求界，保留尺度不确定度下的有效数字。
        at_node = p(arb_poly([x,1]))((scale_squared-1)*x)**2/q(scale_squared*x)
        assert at_node > 0
        lower = (expectation(minorant)/at_node).lower()
        upper = (expectation(squared)/at_node).upper()
        assert 0 < lower < upper
        if i==0:
            assert lower > arb(".21663040965174") and upper < arb(".21663040965176")
        if i==1:
            assert lower > arb(".09793647200697") and upper < arb(".09793647200699")
        intervals.append(lower.union(upper))
        print("零点编号与权重两端",i+1,lower.str(18),upper.str(18))

    cap = (200*scale_squared/acb.zeta_zero(200).imag**2).upper()
    upper = local["weight_cap"](times,cap,41,[0,cap],error)
    assert upper < arb(".121744")
    print("五荷三模总权重上界",upper.str(13))
    return intervals


def large_scale_weights(times, heights, cutoff, error):
    moments = clock["cleared_moments"](times,40)
    gram = arb_mat([[moments[i+j] for j in range(20)] for i in range(20)])
    coefficients = gram.solve(arb_mat([[-moments[20+i]] for i in range(20)]))
    # 首一正交多项式的二十个根都在目标测度支撑的凸包内。
    orthogonal = arb_poly([coefficients[i,0] for i in range(20)]+[1])
    numerator = orthogonal**2
    poles = local["partial_fractions"](40)
    q = arb_poly([1])
    for rate,_ in poles:
        q *= arb_poly([1,arb(rate)/200])
    residues = local["residues_of"](numerator,poles)
    nominal = sum(numerator[i]*moments[i] for i in range(41))
    expectation = nominal+error*sum(abs(c) for c in residues)
    for count,bound in ((1,"7.28e-51"),(199,"1.404e-50"),(8192,"9.07e-25")):
        z = 200*(cutoff/heights[count-1])**2
        assert z > 200/heights[0]**2 and numerator(z)>0
        upper = expectation*q(z)/numerator(z)
        assert 0 < upper < arb(bound)
        print("大尺度的前若干枚总权重上界",count,upper.str(12))


def verify():
    ctx.prec,ctx.cap = 112,3
    heights = [z.imag for z in acb.zeta_zeros(1,15000)]
    assert len(heights)==15000
    # 单独重验三模首次出现于编号 200，以便使用其高度截断。
    local["mask"]["verify"](200)
    ctx.prec = 1536
    times = clock["times"]()
    cutoff = arb("588.4188")
    error = arb(2)**-176
    windows = local["required_windows"](times)
    scale_squared = identify_pair(heights,windows,cutoff)
    small_scale_weights(times,scale_squared,error)
    large_scale_weights(times,heights,cutoff,error)


if __name__ == "__main__":
    verify()
