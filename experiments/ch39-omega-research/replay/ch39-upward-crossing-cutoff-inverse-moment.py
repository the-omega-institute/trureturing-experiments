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

    # 对不确定质量，令负系数先退出、正系数后退出，得到质量尾和的逐点上界。
    jumps = []
    for b, z in heavy:
        assert b > 0 or b < 0
        jumps.append((z.upper() if b > 0 else z.lower(), b/(4*pi)))
    prefix = arb(0)
    for _, weight in sorted(jumps)[:-1]:
        prefix += weight
        assert prefix > 0
    # 总系数精确为零，由 direction 的第三个定义方程保证。
    slope = -sum(b*z for b, z in heavy)/(4*pi)
    target_log = (arb("1e16")/model["mz"]).log()
    gap = target_log-constant
    assert slope > 0 and gap/slope > arb("12755720.35")

    native = run_path(str(root/"ch39-native-rh-ray-pvm.py"))
    weights = {}
    for charge, _, direct, _ in native["shell_weights"](root):
        if charge in (-4, 2):
            pivot = 0 if charge == -4 else 2
            weights[charge] = [a.real for a in direct[pivot]]
    positive = sum((a*b).max(0) for a in weights[-4] for b in weights[2])
    negative = positive-1

    def response(y):
        # -Σ b Ein(zy)/(4π)，因 Σb=0，公共的 log y 精确抵消。
        return dot(direction, initial)-constant-sum(
            b*lib["expint"](1, z*y) for b, z in heavy
        )/(4*pi)

    # 正部用凹函数的 Jensen 上界，负部用 y=1 处的下界。
    low, high = arb(12782662), arb(12782664)
    cost = negative*response(arb(1))
    assert positive*response(low/positive)-cost < gap
    assert positive*response(high/positive)-cost > gap
    print("固定交汇能标的正部倒数矩下界、指定361权的正质量及凹性门槛",
          (gap/slope).str(15), positive.str(15), low, high)


if __name__ == "__main__":
    verify()
