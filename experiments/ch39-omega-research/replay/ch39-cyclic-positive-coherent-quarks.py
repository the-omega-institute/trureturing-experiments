"""任意三循环不变的四站正算子压缩到共同三源空间时的夸克数据障碍。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def window(center, halfwidth):
    return arb(center)+arb(0,arb(halfwidth).upper())


def extension_angle(gram):
    normal = [arb(1)]
    for i in range(1,4):
        normal.append(-sum(gram[i][j]*normal[j] for j in range(i))/gram[i][i])
    return ((normal[0]**2+sum(normal[1:])**2/3)/sum(v*v for v in normal)).sqrt()


def plane_case(exponent, a, mixing):
    # 三族的 Yukawa 边际区间。两种读法分别取一次和二次幂。
    light,middle,heavy = (window(".00000615",".00000014"),window(".00311",".00005"),window(".8616",".0043"))
    u,s,t = [x**exponent for x in (light,middle,heavy)]
    d,m,b = [x**exponent for x in (window(".0000135",".0000002"),window(".000268",".000003"),window(".01401",".00009"))]
    sectors = [(u,s,t),(d,m,b)]
    radii, hierarchies = [], []
    for first,second,third in sectors:
        ratio = ((first+second)/third).upper()
        hierarchies.append(ratio)
        radii.append((ratio/(1+ratio)).sqrt().asin().upper())
    # 不同类型支撑方向之间的最小夹角。
    fixed_to_fixed = ((1-a*a)/(1+a*a)).acos()
    fixed_to_plane = (a/(1+a*a).sqrt()).acos()
    gap = min(fixed_to_fixed.lower(),fixed_to_plane.lower())
    vtb = abs(mixing[2][2])
    angle = vtb.acos()
    assert sum(radii)+angle < gap
    # 两个重方向若都靠近同一固定射线，则用本征方程进一步收紧近邻角。
    fixed_radii = list(radii)
    for i,(ratio,radius) in enumerate(zip(hierarchies,radii)):
        for _ in range(3):
            upper = (ratio/(gap-radius).sin()**2).asin().upper()
            radius = min(radius,upper)
        fixed_radii[i] = radius
    assert sum(fixed_radii) < angle
    # 两边重方向现在都靠近同一平面。用方差条件缩小法向对重本征线的重叠。
    def heavy_cap(second_third, radius):
        p = radius.sin().upper()
        denominator = 1-second_third-p*p-a*p
        assert denominator > 0
        return (a*second_third/denominator).upper()
    pt,pb = [heavy_cap((second/third).upper(),radius) for (_,second,third),radius in zip(sectors,radii)]
    vtd,vts,vtb = (abs(v) for v in mixing[2])
    cd,cs,cb = (abs(v) for v in mixing[1])
    ps = arb(1)
    for _ in range(4):
        ps = min(ps,((vtd+vtb*pb+pt)/vts).upper())
        pc = min(arb(1),(cd+cs*ps+cb*pb).upper())
        # 下列三式用每一根本征线的 p_i |lambda_i-mean| <= a mean。
        mean_up = (u/t + s/t*pc**2 + pt**2).upper()
        mean_down = (d/b + m/b*ps**2 + pb**2).upper()
        assert mean_up<1 and mean_down<1
        pt = min(pt,(a*mean_up/(1-mean_up)).upper())
        pb = min(pb,(a*mean_down/(1-mean_down)).upper())
        mean_strange = (d/m + ps**2 + b/m*pb**2).upper()
        if mean_strange<1:
            ps = min(ps,(a*mean_strange/(1-mean_strange)).upper())
    assert 1-ps**2-pb**2 > 0
    pd = (1-ps**2-pb**2).sqrt().lower()
    required = ((vtd*pd-vtb*pb-pt)/vts).lower()
    print("目标幂、θtb、固定射线半径和、顶底投影上界、奇投影必要下界与上界",exponent,
          angle.lower().str(9),angle.upper().str(9),sum(fixed_radii).str(9),
          pt.str(9),pb.str(9),required.str(9),ps.str(9),flush=True)
    assert required>ps
    return required,ps


def verify(root):
    native=run_path(str(root/"ch39-native-rh-ray-pvm.py"))
    ckm=run_path(str(root/"ch39-shell-clock-ckm-blocks.py"))
    ctx.prec=320
    a=extension_angle(native["gram_cholesky"](root))
    assert arb(".51")<a<arb(".53")
    print("四站法向的平凡分量平方范数",(a*a).str(14),flush=True)
    book = ckm["target_box"]()
    literature = ckm["standard_matrix"](
        window(".22704",".00082").sin(),window(".04275",".00042").sin(),
        window(".00377",".00008").sin(),window("1.139",".023"))
    for name,matrix,limits in (("书库",book,((".20659",".001395"),(".19072",".028926"))),
                               ("1 TeV 同表",literature,((".19959",".001395"),(".18417",".028926")))):
        print(name)
        for exponent,(lo,hi) in zip((2,1),limits,strict=True):
            required,allowed=plane_case(exponent,a,matrix)
            assert required>arb(lo) and allowed<arb(hi)

    # 正扩张条件只依赖法向在平凡表示中的长度，故可一次覆盖四种选源。
    triplets = run_path(str(root/"ch39-four-gram-triplet-ckm.py"))
    gram = native["gram_cholesky"](root)
    band = arb(".5198").union(arb(".5774"))
    for omit in range(4):
        n = [v.real for v in triplets["normal"](gram,omit)]
        length = (n[0]**2+sum(n[1:])**2/3).sqrt()
        assert band.contains(length)
    for matrix in (book, literature):
        for exponent in (1,2):
            required,allowed = plane_case(exponent,band,matrix)
            assert required > arb(".18") and allowed < arb(".034")


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
