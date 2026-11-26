"""105圈固定真空沿味读数方向的局部最小能量代价，与实际零点的比较。"""
from pathlib import Path
from runpy import run_path

from flint import arb, acb, ctx


def verify():
    root = Path(__file__).resolve().parent
    fixed = run_path(str(root/"ch39-height-independent-flavor-vacuum.py"))
    angles = run_path(str(root/"ch39-mincut-reference-ckm-action.py"))
    winding = run_path(str(root/"ch39-mincut-global-winding.py"))
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    cp = run_path(str(root/"ch39-cp-ensemble-test.py"))
    ctx.prec, ctx.cap = 384, 3
    flavor = cp["flavor"]
    d = angles["reference_angles"](flavor)[0]
    x, y = arb("-.01046"), arb("-.01046")-d
    tri = winding["terms"](shell)
    weights = fixed["minimum_order_weights"](winding, tri, x, y)
    cap = arb(23)/25*arb(2)**-20/14**2
    modes = [(q, [w*z for z in winding["coeff"](tri, q, arb(0).union(cap))])
             for q, w in weights.items()]
    _, _, a, b, h = winding["derivatives"](modes, x, y)
    # 最小化 v^T H v，约束 v_x-v_y=1。
    k = (a*h-b*b)/(a+2*b+h)
    ax, bx = [flavor["singlet_at"]({724:z,804:acb(1)})[1] for z in (acb(1),acb(0,1))]
    c = flavor["singlet"](11661)[1]
    data = cp["two"]["spectral_data"]()
    polys = cp["cp_data"](data, flavor["reference_matrix"]())
    u = cp["two"]["probabilities"](c, data)[0]
    j = cp["cp_value"](c, polys, True)
    slope = -ax*d.sin()+bx*d.cos()
    stiffness = [k/(slope*v[1]/v[0])**2 for v in (u,j)]
    print("钟差、相对 U 和相对 J 的局部最优二次系数（未除以二）",
          *[v.str(16) for v in (k,*stiffness)])
    assert arb("2229.007") < k < arb("2229.009")
    assert arb("23.40204") < stiffness[0] < arb("23.40206")
    assert arb("62.76770") < stiffness[1] < arb("62.76774")
    homogeneous = (a+2*b+h)/2
    assert arb("1026.3462") < homogeneous < arb("1026.3464")
    print("同时等角推动两钟的二次能量系数", homogeneous.str(16))


if __name__ == "__main__":
    verify()
