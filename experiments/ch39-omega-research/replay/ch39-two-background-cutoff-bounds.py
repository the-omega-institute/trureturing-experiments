"""任意归一有符号截止下，两个背景共用真空读数及电磁—色匹配的必要条件。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify():
    root = Path(__file__).resolve().parent
    fields = run_path(str(root / "ch39-two-zero-cutoff-gravity.py"))["gravity_fields"]
    vector_moment = run_path(str(root / "ch39-em-colour.py"))["vector_log_moment"]
    _, _, (rows,) = fields(root, (2_306_888_100_085,))
    masses = [z for _, _, z in rows[49:65]]
    assert len(masses) == 16
    ctx.prec = 256
    v, planck, mz = map(arb, ("246.22", "1.22089e19", "91.1876"))
    g1, g2, yt = map(arb, (".467453", ".63811", ".8616"))
    pi = arb.pi()
    r = (v/planck)**2/2
    x_min = arb("4e-15")
    k_w, k_z, k_t = g2**2/2, (g2**2+3*g1**2/5)/2, yt**2
    assert 0 < k_w < k_z < k_t and min(masses) > r
    inverse_masses = sum(1/(z-r) for z in masses)
    scalar_charge = r*inverse_masses/(12*pi)
    heavy_drop = r*inverse_masses/(48*pi**2)
    flat = arb(".129")
    light =  ((12, k_t), (-6, k_w), (-3, k_z))
    max_light = sum(n*k**2 for n,k in light)/(128*pi**2)
    light_loss = (r/x_min)*sum(n*k**3 for n,k in light)/(192*pi**2)
    unrestricted_tv = (flat+max_light/2)/(max_light/2+heavy_drop)
    supported_tv = (flat+max_light-light_loss/2)/(light_loss/2+heavy_drop)
    min_drop = 8*pi*r**2*(flat-12*max_light-25*heavy_drop)
    assert unrestricted_tv > arb("54.02266")
    assert supported_tv > arb("2.0766053e21")
    assert min_drop > arb("7.33892e-68")
    print("任意正节点的必要总变差、固定支撑下限的必要总变差、总变差≤25时的真空落差下界",
          unrestricted_tv.lower().str(17), supported_tv.lower().str(17), min_drop.lower().str(17))

    z_w = r*k_w
    p0, min_p = (mz/planck)**2, (arb("1e13")/planck)**2
    factor = 7/(4*pi)
    reference = factor*vector_moment(z_w/p0)
    tail = factor*vector_moment(z_w/min_p)
    unit = factor*z_w/(2*x_min)
    colour_tv = (reference-tail-unit)/(unit+scalar_charge)
    error = tail+26*unit+25*scalar_charge
    assert colour_tv > arb("3.5087950e20")
    assert error < arb("7.76e-20")
    print("电磁与色的位移参考、总变差≤25的误差上界、双背景同匹配的必要总变差",
          (-reference).str(26), error.upper().str(15), colour_tv.lower().str(18))


if __name__ == "__main__":
    verify()
