"""给出实际零点同时逼近105圈真空作用与两个味读数的显式相位窗和下密度。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def verify(root=None):
    root = Path(__file__).resolve().parent if root is None else root
    fixed = run_path(str(root / "ch39-height-independent-flavor-vacuum.py"))
    angles = run_path(str(root / "ch39-mincut-reference-ckm-action.py"))
    winding = run_path(str(root / "ch39-mincut-global-winding.py"))
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    cp = run_path(str(root / "ch39-cp-ensemble-test.py"))
    ctx.prec, ctx.cap = 384, 3
    flavor = cp["flavor"]
    d = angles["reference_angles"](flavor)[0]
    x, y = arb("-.01046"), arb("-.01046")-d
    angle_radius = arb("9e-6")
    difference_radius = arb("2e-7")
    cap = arb(23)/25*arb(2)**-20/14**2
    tri = winding["terms"](shell)
    weights = fixed["minimum_order_weights"](winding, tri, x, y)
    modes = [(q, [w*z for z in winding["coeff"](tri, q, arb(0).union(cap))])
             for q, w in weights.items()]
    # 新坐标 s,t 给原钟位移 (s+t,s)。目标处两个梯度严格为零。
    xx, yy = x+arb(0, angle_radius+difference_radius), y+arb(0, angle_radius)
    _, _, a, b, h = winding["derivatives"](modes, xx, yy)
    energy_cap = (abs(a+2*b+h)*angle_radius**2+2*abs(a+b)*angle_radius*difference_radius
                  + abs(a)*difference_radius**2)/2
    assert energy_cap < arb("1e-7")
    ax, bx = [flavor["singlet_at"]({724:z,804:acb(1)})[1] for z in (acb(1),acb(0,1))]
    t = d+arb(0,difference_radius)
    coordinate = ax*t.cos()+bx*t.sin()
    reference = flavor["singlet"](11661)[1]
    data = cp["two"]["spectral_data"]()
    polys = cp["cp_data"](data, flavor["reference_matrix"]())
    # 在整个窄角差段用导数界，避免概率的整式区间包围损失小量精度。
    c_error = abs(-ax*t.sin()+bx*t.cos())*difference_radius
    u0 = cp["two"]["probabilities"](reference, data)[0][0]
    j0 = cp["cp_value"](reference, polys)
    u_error = abs(cp["two"]["probabilities"](coordinate, data)[0][1])*c_error/u0
    j_error = abs(cp["cp_value"](coordinate, polys, True)[1])*c_error/j0
    assert u_error < arb("2e-6") and j_error < arb("1.2e-6")
    # 谱的变分用本征值范数界，避免对整个角度区间直接用高阶矩阵求根。
    physical = run_path(str(root / "ch39-first-three-mode-observables.py"))
    model, initial = physical["model_inputs"]()
    chiral = physical["joint"]["chiral"]
    read = chiral["source"]["spectral_supplier"](model)
    levels = read(arb("1e13"), clocks={724:acb(0,x).exp(), 804:acb(0,y).exp()})
    shell_cap = ((2*36+15)*4+36/3)*cap
    # 利用三特征标在整个相位方框上给出谱区间。
    inertia = run_path(str(root/"ch39-neutral-mean-phase-inertia.py"))
    spectral = run_path(str(root/"ch39-twentyseven-coefficient-window.py"))
    neutral = chiral["ordered"]([v for v, q in levels if q == 0])
    neutral = [v+arb(0,shell_cap) for v in neutral]
    delta = arb(0).union(cap)
    bands = {q:spectral["bands"](shell, inertia, q, xx, yy, delta) for q in model["charges"]}
    bands[0] = neutral
    assert all(bands[q][2]<9 and bands[q][3]>10 for q in model["charges"])
    _, alpha, curve = physical["layer_observables"](model, initial, bands)
    assert arb(".117")<alpha<arb(".119") and arb(".998")<curve<arb("1.002")
    density = angle_radius*difference_radius/arb.pi()**2
    assert density > arb("1.82378e-13")
    print("共同开窗的相对下密度、能差上界、U 与 J 相对误差上界",
          *[v.str(15) for v in (density, energy_cap, u_error, j_error)])
    return model, initial, physical["joint"], bands, density, delta


if __name__ == "__main__":
    verify()
