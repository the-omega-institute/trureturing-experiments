"""允许任意字符保持门时，实现固定六壳见证至少需要多少次原始 Gram 射线相位门。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def experiment(root):
    def lib(name):
        return run_path(str(root/("ch39-" + name + ".py")))

    gram = lib("native-rh-ray-pvm")["gram_cholesky"](root)
    ctx.prec = 256
    angles = [arb(0)]
    assert all(value.is_zero() for value in gram[0][1:])
    for pole, first, second, third in gram[1:]:
        squared = pole**2 + first**2 + second**2 + third**2
        weight = (pole**2 + (first+second+third)**2/3)/squared
        assert 0 < weight < 1
        mix = 2*(weight*(1-weight)).sqrt()
        assert mix < arb(".18283")
        angles.append(mix.asin())
    largest = max(angle.upper() for angle in angles)
    # 公开主角证书给出固定六壳认证盒内每个根的所有实现酉都超过此量。
    required = arb(".941781").asin()
    assert 6*largest < required
    print("单个原始 Gram 相位门的最大主角上界、见证主角下界（度）",
          (largest*180/arb.pi()).str(12), (required*180/arb.pi()).str(12))


if __name__ == "__main__":
    experiment(Path(__file__).resolve().parent)
