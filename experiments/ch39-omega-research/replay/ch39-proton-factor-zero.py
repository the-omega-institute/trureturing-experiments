"""保持同一跨高度参考味真空时，十二面正反等分的质子因子可精确为零。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx


# 以下十进有理数来自固定真空与质子因子的联合选权。
FIXED = {
    1: ".017671857299", 3: ".020134166880", 5: ".001158304825",
    10: ".001931646563", 12: ".001239021195", 14: ".001749647913",
    16: ".000350538208", 32: ".022095550408", 63: ".000630825890",
    65: ".001369477075", 67: ".027618176976", 69: ".008222871705",
}
SOLVED = (33, 35, 36, 38, 104, 105)


def verify():
    root = Path(__file__).resolve().parent
    fixed = run_path(str(root / "ch39-height-independent-flavor-vacuum.py"))
    global_model = run_path(str(root / "ch39-reference-flavor-global-winding.py"))
    lib = run_path(str(root / "ch39-mincut-global-winding.py"))
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    angles = run_path(str(root / "ch39-mincut-reference-ckm-action.py"))
    flavor = run_path(str(root / "ch39-chiral-neutral-flavor-phase.py"))
    ctx.prec = 224
    x = arb("-.01046")
    difference = angles["reference_angles"](flavor)[0]
    y = x - difference
    tri = lib["terms"](shell)
    global_model["conjugation_symmetry"](root)
    assert len(tri[1, -1]) == 12 and all(n == 0 for _, n in tri[1, -1])
    height = acb.zeta_zero(2_306_888_100_085).imag
    phi = shell["geo"]["phases"](height)
    actual = (phi[724]/phi[804]).arg()
    cap = arb(23)/25*arb(2)**-20/14**2
    delta = arb(0).union(cap)
    solutions = []
    for cancellation, other in ((difference, actual), (actual, difference)):
        weights = {q: arb(v) for q, v in FIXED.items()}
        bare = {q: lib["coeff"](tri, q, arb(0))[2]/12 for q in (*FIXED, *SOLVED)}

        def factor(q, t):
            return (bare[q]*acb(0, q*t).exp()).real

        columns = {q: [arb(1), *fixed["gradient_columns"](lib, tri, q, x, y),
                       factor(q, cancellation)] for q in bare}
        system = arb_mat([columns[q] for q in SOLVED]).transpose()
        rhs = [int(j == 0) - sum(w*columns[q][j] for q, w in weights.items())
               for j in range(6)]
        assert not system.det().contains(0)
        values = system.solve(arb_mat([[b] for b in rhs])).entries()
        weights.update(zip(SOLVED, values, strict=True))
        assert min(weights.values()) > arb(".00035")
        modes = [(q, [w*c for c in lib["coeff"](tri, q, delta)]) for q, w in weights.items()]
        drift = global_model["global_cover"](lib, tri, weights, modes, x, y)
        assert drift < -3
        value = sum(w*factor(q, other) for q, w in weights.items())
        gradient = sum(w*(acb(0, q)*bare[q]*acb(0, q*cancellation).exp()).real
                       for q, w in weights.items())
        window = difference + arb(0, arb("2e-7"))
        derivative = sum(w*(acb(0, q)*bare[q]*acb(0, q*window).exp()).real
                         for q, w in weights.items())
        # 同编号自身相位与同高度自由真空的能量差，使用本次解得的权。
        dnumber = arb(23)/25*arb(2)**-20/(height**2+shell["soft"]["point"]["study"]["offset"])
        energy = arb(0)
        for q, w in weights.items():
            cs = lib["coeff"](tri, q, dnumber)
            energy += w*(lib["energy"](q, cs, phi[724].arg(), phi[804].arg())
                         - lib["energy"](q, cs, x, y))
        print("相消角差、另一处读数、相消点导数、整窗导数、实际零点激发能",
              *(v.str(16) for v in (cancellation, value, gradient, derivative, energy)))
        assert not value.contains(0) and abs(value) < arb("1e-6")
        assert 0 < energy < arb("1e-7")
        assert arb("1.069") < derivative < arb("1.071")
        assert sum(w*factor(q, cancellation) for q, w in weights.items()).contains(0)
        solutions.append((value, energy))
        if cancellation is difference:
            assert value < 0
        else:
            assert value > 0


    # 一组新权仍跨高度固定：把上面两个物理解沿直线混合，
    # 令两个指定角差上的因子互为相反数。激发能也沿同一直线变化。
    numbered_in_first, first_energy = solutions[0]
    vacuum_in_second, second_energy = solutions[1]
    share = vacuum_in_second/(vacuum_in_second-numbered_in_first)
    residual = (1-share)*vacuum_in_second
    mixed_energy = share*first_energy+(1-share)*second_energy
    assert 0 < share < 1 and 0 < residual < arb("7.28e-8")
    assert 0 < mixed_energy < arb("1e-7")
    print("使两点残差等模反号的第一组权份额与最大模", share.str(17), residual.str(16))

if __name__ == "__main__":
    verify()
