"""上行一圈交汇对归一有符号截止倒数矩的必要条件。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx


def verify():
    root = Path(__file__).resolve().parent
    fields = run_path(str(root/"ch39-two-zero-cutoff-gravity.py"))["gravity_fields"]
    _, lib, (central,) = fields(root, (2_345_666,), tight_neutral=True)
    model, _, _, initial, _ = lib["chiral"]["source"]["mean"]["inputs"]()
    ctx.prec = 192
    pi = arb.pi()
    beta = [v/(2*pi) for v in model["beta"]]
    coefficients = run_path(str(root/"ch39-signed-gravity-unification.py"))["gauge_weights"](model)
    total = [sum(b[i] for b in coefficients) for i in range(3)]
    equations = arb_mat([[1, 1, 1], beta, total])
    assert not equations.det().contains(0)
    direction = equations.solve(arb_mat([[0], [1], [0]])).entries()
    dot = run_path(str(root/"ch39-five-zero-gravity-gauge-transfer.py"))["dot"]
    heavy = [(dot(direction, b), z) for b, (_, _, z) in zip(coefficients, central[1:], strict=True)]
    constant = dot(direction, initial)+sum(b*(arb.const_euler()+z.log()) for b, z in heavy)/(4*pi)
    envelope = sum(abs(b)*z for b, z in heavy)/(4*pi)
    floor = -constant/envelope
    print("消去共同逆耦合和共同对数后的常数、余项系数、上行交汇的截止倒数矩下界",
          constant.str(13), envelope.str(13), floor.str(13))
    assert constant < -8344 and envelope > 0 and floor > arb("12671470.84")


if __name__ == "__main__":
    verify()
