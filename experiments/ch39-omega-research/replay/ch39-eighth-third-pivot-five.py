"""第八壳实际第三枢轴的五个荷六端点条件与原源谱主序。"""
from pathlib import Path
from runpy import run_path


def verify(root):
    solve = run_path(str(root / "ch39-fifth-third-pivot-em-colour-path.py"))["solve"]
    centers = (
        "0.004932616450557562",
        "7.888943691425847e-16",
        "2.6468714117795603e-11",
        "9.843149807488264e-17",
        "3.2420032403979707e-25",
        "0.0002556936313929033",
        "1.8084969303171906e-18",
        "0.24984786294117972",
        "7.284750387480544e-12",
        "1.137242633441764e-34",
        "0.05786442774855582",
        "0.09604716399622006",
        "2.4328998424440664e-21",
        "1.6343601089703911e-15",
        "2.604341804006931e-24",
        "3.887828103050434e-10",
        "0.0023538334475496823",
        "3.753952530442964e-11",
        "5.292989308042476e-12",
    )
    data, rows, weights, sites, boxes, _, lib = solve(root, centers, (9, 16, 18, 7, 15, 0), number=2_306_888_100_093)
    from flint import arb
    gauge = lib("all-field-signed-crossing")
    full = lib("full-neutral-running")
    _, _, color = full["direct"](data, sites, weights, rows, "560000", gauge, gauge["angular_rule"](16,45))
    norm = sum(abs(arb(coeff[2])) for _, _, coeff in rows)/(4*arb.pi())
    error = norm*sum(abs(w)*abs((b/x).log()) for w,b,x in zip(weights,boxes,sites,strict=True))
    print("该精确根的 56 万 GeV 色逆耦合上界", (color+error).upper().str(13),flush=True)
    assert color+error < -30



if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
