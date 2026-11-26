"""三射线相对重相的整个相位环面与四种 Gram 三源的正夸克读出。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def geometry_floor(normal):
    ray_moduli = [abs(z) for z in normal[1:]]
    floor = (2*max(ray_moduli)-sum(ray_moduli))/arb(3).sqrt()
    assert 0 < floor < 1
    # 三个 Fourier 分量模长的反三角界，加上 a²+b²+c²=1。
    slope_cap = (1-floor**2)/(2*floor)
    pair_gap = 2*floor.atan()
    plane_gap = (1/slope_cap).atan()
    return floor,slope_cap,min(pair_gap.lower(),plane_gap.lower())


def verify(root):
    native = run_path(str(root/"ch39-native-rh-ray-pvm.py"))
    triplets = run_path(str(root/"ch39-four-gram-triplet-ckm.py"))
    real_cone = run_path(str(root/"ch39-cyclic-positive-coherent-quarks.py"))
    ckm = run_path(str(root/"ch39-shell-clock-ckm-blocks.py"))
    ctx.prec = 320
    gram = native["gram_cholesky"](root)
    slope = arb(".73")
    # 旧证书只用方差系数上界和不同生成成分的角距下界。
    assumed_gap = min((2*slope.atan()).lower(),(1/slope).atan().lower())
    geometry = ((".5084305",".729204",".9407385"),(".5339539",".669434",".9808807"),
                (".5570340",".619095","1.0164551"),(".5773502",".577351","1.0471975"))
    for omitted, (low,high,separation) in enumerate(geometry):
        normal = triplets["normal"](gram,omitted)
        floor,cap,gap = geometry_floor(normal)
        assert floor > arb(low) and cap < arb(high) < slope
        assert gap > arb(separation) > assumed_gap
        print("省略的 Gram 源、三类分量模长下界、斜率上界与角距下界",
              omitted,floor.lower().str(13),cap.upper().str(13),gap.str(13))

    win = real_cone["window"]
    same_scale = ckm["standard_matrix"](
        win(".22704",".00082").sin(),win(".04275",".00042").sin(),
        win(".00377",".00008").sin(),win("1.139",".023"))
    for name,matrix,limits in (("书库",ckm["target_box"](),((".20658",".001961"),(".17917",".047039"))),
                              ("1 TeV 同表",same_scale,((".19958",".001961"),(".17278",".047306")))):
        print(name)
        for exponent,(lo,hi) in zip((2,1),limits,strict=True):
            required,allowed = real_cone["plane_case"](exponent,slope,matrix)
            assert required>arb(lo) and allowed<arb(hi)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
