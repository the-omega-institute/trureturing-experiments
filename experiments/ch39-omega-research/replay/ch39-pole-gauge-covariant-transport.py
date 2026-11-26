"""用原壳的十八条极点辐边确定跨壳平凡化，并核对两端独立换规。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def adjoint(a):
    return a.conjugate().transpose()


def diagonal(values):
    return acb_mat([[s if i == j else 0 for j in range(len(values))]
                    for i, s in enumerate(values)])


def spoke_frame(laplacian, pole):
    return diagonal([1 if i == pole else -laplacian[i, pole] / abs(laplacian[i, pole])
                     for i in range(19)])


def verify(root):
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    pvm = run_path(str(root / "ch39-adjacent-pvm-transport.py"))
    ctx.prec = 512
    pole = shell["soft"]["pole"]
    assert pole == 14
    center = 2_306_888_100_085
    identity = diagonal([1]*19)
    for q in (-4, 2, 6):
        prior = None
        for offset in (-1, -2, -3, 1, 2, 3):
            n = center + offset
            h = shell["laplacian"](q, acb.zeta_zero(n).imag)
            assert all(h[i, pole].real < 0 and h[i, pole].imag.is_zero()
                       for i in range(19) if i != pole)
            assert spoke_frame(h, pole).contains(identity)
            _, basis = pvm["spectrum"](shell, q, n)
            # 仅作复算示例。证明适用于两端各站任意独立相位。
            phase = []
            for i in range(19):
                k = acb(0, arb(i + 1)*offset/19)
                phase.append(((1+k)/(1-k))**q)
            gauge = diagonal(phase)
            corrected_frame = spoke_frame(gauge*h*adjoint(gauge), pole)
            corrected_basis = adjoint(corrected_frame)*gauge*basis
            if abs(offset) > 1:
                before, after = prior
                reference = adjoint(basis)*before
                measured = adjoint(corrected_basis)*after
                errors = [abs(a)**2 - abs(b)**2 for a, b in
                          zip(reference.entries(), measured.entries(), strict=True)]
                assert all(d.contains(0) for d in errors)
                width = max(abs(d).upper() for d in errors)
                print(q, offset, width.str(5))
                assert width < arb("1e-35")
            prior = basis, corrected_basis
    print("三荷各四步、每步 361 个投影重叠在两端独立换规后核验通过")


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
