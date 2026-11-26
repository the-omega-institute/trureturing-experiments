"""在全双钟域给同一非负固有时轮廓的重规范阈值及三支逆耦合残差作界。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify():
    root = Path(__file__).resolve().parent
    gravity = run_path(str(root / "ch39-chiral-induced-gravity-vacuum.py"))
    spectrum = gravity["spectrum"]
    model, _, _, initial, _ = spectrum["torus"]["local"]["inputs"]()
    ctx.prec = 224
    bands = spectrum["shell_bands"](model)
    assert gravity["prefix_cap"](bands, 32) == 0
    assert all(h > 0 for values in bands.values() for h in values)
    fields = run_path(str(root / "ch39-chiral-minimal-neutral-gravity.py"))["fields"]
    lower = {q: [h.lower() for h in values] for q, values in bands.items()}
    rows = fields(lower)
    weights = []
    for q, coeff in zip(model["charges"], model["weyl"], strict=True):
        beta = [2*f + (model["higgs"][a] if q == 3 else 0) for a, f in enumerate(coeff)]
        weights.extend([beta]*16)
    weights.extend([[0, -7, -arb(21)/2]]*18)
    assert len(weights) == len(rows)-1 == 98
    epsilon = arb("1e-122")
    limits = [sum(abs(b[a]) / z.lower()**2 for b, (_, _, z) in zip(weights, rows[1:], strict=True))
              * 2*epsilon/(31*arb(2).exp()) for a in range(3)]
    assert all(v < arb("1e-99") for v in limits)
    initial[2] = 1/arb(".117").union(arb(".119"))
    b1,b2,b3 = model["beta"]
    direction = [b2-b3, b3-b1, b1-b2]
    norm = sum(abs(u) for u in direction)
    # 该方向与共同匹配值及一圈对数斜率均正交。
    spread = 2*(abs(sum(d*a for d, a in zip(direction, initial)))
                -sum(abs(d)*v for d, v in zip(direction, limits)))/norm
    assert spread > arb("3.61")
    print("全双钟域三支重阈值绝对值上界",*[v.upper().str(12) for v in limits])
    print("任意能标的最大逆耦合减最小逆耦合的下界",spread.lower().str(15))


if __name__ == "__main__":
    verify()
