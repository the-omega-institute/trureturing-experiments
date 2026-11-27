"""第八壳弧上正色根在第九壳实际谱基中的两个引力读数。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify(root):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    arc = lib("eighth-shortest-arc-pivot-support")
    eighth = arc["NUMBER"]
    ninth = eighth + 1
    boxes = lib("eighth-arc-positive-probe")["verify"](root)
    rows, single = arc["polynomials"](root, number=ninth)
    weights = (single[1], arc["angle_weights"](rows[1], arb(3) / 4))
    shell = lib("three-ray-chiral-boundary")
    spectrum = lib("adjacent-pvm-transport")["spectrum"]
    u, v = [spectrum(shell, 6, n)[1] for n in (eighth, ninth)]
    overlap = v.conjugate().transpose() * u
    moments, heat, (field,) = lib("three-index-free-atoms")["broken_fields"](root, [ninth])
    _, source = lib("adjacent-ordered-charge-six-cutoffs")["verify"](
        root, refine_node=2_306_888_100_086
    )[2]
    ctx.prec = 192
    cap = arb(".75553")
    assert max(source).sqrt() < cap
    slope_k = sum(a for _, a, _ in field) / (2 * arb.pi())
    slope_w = -sum(g for g, _, _ in field) / (8 * arb.pi())
    errors = lib("fifth-seventh-cutoff-separation")["moment_errors"]
    for case, (box, row) in enumerate(zip(boxes, weights, strict=True)):
        # 所有在第八壳具有这组对角元的正算子均受下列界约束。
        bounds = []
        for j in range(19):
            terms = [abs(overlap[j, i]) * y.sqrt() for i, y in enumerate(box)]
            upper = sum(terms).upper()**2
            lower = max((t - sum(terms[:i] + terms[i + 1:])).lower() for i, t in enumerate(terms))
            bounds.append((lower.max(arb(0))**2).union(upper.min(cap)))
        mk, mw = errors(field, heat, row, [".75553"] * 19)
        normal = [sum(w * ys.nonnegative_part()**p for w, ys in zip(row, bounds, strict=True))
                  for p in (2, 4)]
        k = slope_k * (normal[0] + (mk - 1 / slope_k))
        w = slope_w * (normal[1] + (mw - arb("1e-122") / slope_w))
        ordered = [moments(field, heat, y**2) for y in box]
        by_index = [sum(c * m[j] for c, m in zip(row, ordered, strict=True)) for j in (0, 1)]
        print("读出编号", case, "同编号复用的K、W", *[z.str(10) for z in by_index],
              "同算子的K、W外包", *[(z.lower().str(10), z.upper().str(10)) for z in (k, w)], flush=True)
        k_lo, k_hi, w_lo, w_hi = (
            ("-.08581", "-.07116", ".10838", ".11088"),
            (".02653", ".07374", "-.04731", "-.01071"),
        )[case]
        assert arb(k_lo) < k < arb(k_hi) and arb(w_lo) < w < arb(w_hi)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
