"""三条软辐边对二十七圈真空的符号选择及原味模型的留出读数。"""
from collections import Counter
from pathlib import Path
from runpy import run_path
from flint import acb, arb, ctx


def verify():
    root = Path(__file__).resolve().parent
    faces = run_path(str(root / "ch39-mincut-real-triangle-tomography.py"))
    global_model = run_path(str(root / "ch39-mincut-global-winding.py"))
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    flavor_model = run_path(str(root / "ch39-three-window-two-probability-cp.py"))["cp"]
    ctx.prec, ctx.cap = 384, 3
    table = faces["enumerate_faces"]()
    tri = global_model["terms"](shell)
    radius = arb(0, "2.1e-17")
    x, y = [arb(mid)+radius for mid in (".0707075643615173", ".07570197643661009")]
    split = arb(0)
    for mode, angle in (((1,0),x), ((0,1),y), ((1,-1),x-y)):
        # 整数对 (a,b) 表示相位方向 a+i b sqrt(867)。严格配对给出单位辐边时的精确镜面对称。
        counts = Counter(z for k,z in table.values() if k==mode)
        assert counts == Counter({faces["conjugate"](z): n for z,n in counts.items()})
        soft = [z for z, n in tri[mode] if n == 1]
        for power, weight in ((1,100), (27,13)):
            imaginary = sum((z**power).imag for z in soft)
            split -= (arb(2)*weight/113)*imaginary*(power*angle).sin()
    assert arb("-.052468845665") < split < arb("-.052468845663")
    print("正角能量减镜像能量除以辐边亏额", split.str(17))

    f = flavor_model["flavor"]
    data = flavor_model["two"]["spectral_data"]()
    polynomials = flavor_model["cp_data"](data, f["reference_matrix"]())
    def readings(c):
        u, v = flavor_model["two"]["probabilities"](c, data)
        return u[0], v[0], flavor_model["cp_value"](c, polynomials)

    coordinate = f["singlet_at"]({724:acb(0,x).exp(), 804:acb(0,y).exp()})[1]
    reference = readings(f["singlet"](11661)[1])
    ratios = [a/b for a,b in zip(readings(coordinate),reference,strict=True)]
    limits = ((".26916424967",".26916424970"), ("1.00498810852","1.00498810853"),
              (".31352152717",".31352152720"))
    assert all(arb(lo)<v<arb(hi) for v,(lo,hi) in zip(ratios,limits,strict=True))
    print("六面味坐标与 U、C、J 相对参考的比值", coordinate.str(16), *[v.str(16) for v in ratios])


if __name__ == "__main__":
    verify()
