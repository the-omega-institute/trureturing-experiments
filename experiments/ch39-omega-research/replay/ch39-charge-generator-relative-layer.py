"""荷六匹配约束同一生成元的两种夸克输运在荷六目的第十六态上的夹角。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


def verify(root):
    train = run_path(str(root / "ch39-three-shell-two-space-readout.py"))["verify"]
    records, _ = train(root, return_data=True)
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    transport = run_path(str(root / "ch39-adjacent-pvm-transport.py"))
    center = 2_306_888_100_085
    ctx.prec = 320
    # 两个概率上界来自 ch39-charge-sixteenth-spectral-window.py 的 K=1、D=0 证书。
    for side, high_count, cap in ((-1, 4, "1.2574e-17"), (1, 2, "1.0131e-17")):
        energy, source = transport["spectrum"](shell, 6, center + side)
        _, dest = transport["spectrum"](shell, 6, center + 2*side)
        weights = transport["training_values"](energy, records)
        high = [j for j, value in enumerate(weights) if value > arb(".001")]
        low = [j for j, value in enumerate(weights) if value < arb(".001")]
        assert len(high) == high_count and len(low) + len(high) == 19
        assert all(weights[j].overlaps(weights[low[0]]) for j in low)
        v = acb_mat([[dest[j, 15]] for j in range(19)])
        coords = source.conjugate().transpose()*v
        b = sum(abs(coords[j, 0])**2 for j in high)
        eps = arb(cap).sqrt()
        # b 是 v 在源高层上的概率，精确解给出的 w=U*v 在源高层上至多有 eps 的范数。
        angle_floor = b.sqrt().asin() - eps.asin()
        loss_floor = angle_floor.sin()**2
        # 在源等幅子空间中对齐 v,w。低空间一项已给出重叠的下界。
        retained = ((1-b)*(1-eps**2)).sqrt()
        assert retained > b.sqrt()*eps
        loss_cap = 1-retained**2
        assert 0 < loss_floor < loss_cap < 1 and loss_cap-loss_floor < arb("3e-9")
        print("源偏移、最小两荷失配概率的下界和上界、最小夹角的度数包围",
              side, loss_floor.str(18), loss_cap.str(18),
              (angle_floor.union(retained.acos())*180/arb.pi()).str(13), flush=True)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
