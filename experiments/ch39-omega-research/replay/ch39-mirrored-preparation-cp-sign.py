"""两枚编号置换与两钟独立反号时，计算正交三源 Gram 对的三次交换子迹。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb


def gram(model, charge, height, orientation):
    adj, shell = model["adj"], model["shell"]
    original = shell["geo"]["phases"](height)
    assert sorted(original) == [724, 804]
    phases = {distance: original[distance] if sign == 1 else original[distance].conjugate()
              for distance, sign in zip((724, 804), orientation, strict=True)}
    laplacian = shell["laplacian"](charge, height, phases=phases)
    rows = []
    for frame in model["frames"]:
        levels, vectors = model["eig"](adj(frame) * laplacian * frame)
        assert levels[0] < 9 < levels[1]
        first = acb_mat([[vectors[i, 0]] for i in range(vectors.nrows())])
        rows.append((adj(frame*first)*model["sources"]/(1+(levels[0]/9)**512)).entries())
    matrix = acb_mat(rows)
    return adj(matrix)*matrix


def verify():
    model = run_path(str(Path(__file__).with_name("ch39-two-zero-quark-window.py")))
    heights = [acb.zeta_zero(index).imag for index in (2345666, 964890)]
    orientations = ((1, 1), (-1, -1), (-1, 1), (1, -1))
    predictions = ((-1.465056, 1.465056, -40.3106, 40.3106), (-1.25908, 1.25908, -85.0639, 85.0639))
    for swap, expected in enumerate(predictions):
        for signs, center in zip(orientations, expected, strict=True):
            up = gram(model, -4, heights[swap], signs)
            down = gram(model, 2, heights[1-swap], signs)
            comm = up*down - down*up
            chirality = (comm*comm*comm).trace().imag * 10**12
            assert abs(chirality - arb(str(center))) < arb(".001")
            print("交换编号、两钟取向、交换子三次迹虚部乘以万亿", swap, signs, chirality.str(10))


if __name__ == "__main__":
    verify()
