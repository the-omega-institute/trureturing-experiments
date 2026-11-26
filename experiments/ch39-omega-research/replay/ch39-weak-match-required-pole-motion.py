"""任意跨壳酉联络满足第十六层低截止所需的极点主角及转移概率。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify(root):
    train = run_path(str(root / "ch39-three-shell-two-space-readout.py"))["verify"]
    records, _ = train(root, return_data=True)
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    transport = run_path(str(root / "ch39-adjacent-pvm-transport.py"))
    center = 2_306_888_100_085
    pole = shell["soft"]["pole"]
    assert pole == 14
    ctx.prec = 320
    cutoff = 16 * arb("91.1876") / arb("1.22089e19")
    for side in (-1, 1):
        source = center + side
        destination = center + 2*side
        energy, vec = transport["spectrum"](shell, 6, source)
        _, out = transport["spectrum"](shell, 6, destination)
        values = transport["training_values"](energy, records)
        low = [j for j, y in enumerate(values) if y < cutoff]
        high = [j for j, y in enumerate(values) if y > cutoff]
        assert len(low) + len(high) == 19 and len(high) == (4 if side == -1 else 2)
        low_pole = sum(abs(vec[pole, j])**2 for j in low)
        high_pole = sum(abs(vec[pole, j])**2 for j in high)
        assert (low_pole + high_pole).contains(1)
        source_floor = min(y.lower() for y in values)
        high_floor = min(values[j].lower() for j in high)
        max_high = ((cutoff-source_floor)/(high_floor-source_floor)).sqrt()
        overlap_upper = low_pole.sqrt() + high_pole.sqrt()*max_high
        target_angle = abs(out[pole, 15]).acos()
        angle_lower = overlap_upper.acos() - target_angle
        angle_upper = low_pole.sqrt().acos() - target_angle
        assert 0 < angle_lower < angle_upper < arb.pi()/2
        assert angle_upper - angle_lower < arb("1e-6")
        print(source, destination, "高层及极点在高层的概率", [j+1 for j in high],
              high_pole.str(12), "低截止所需最小极点转角(度)与站外概率",
              (angle_lower.union(angle_upper)*180/arb.pi()).str(12),
              (angle_lower.sin().union(angle_upper.sin())**2).str(12), flush=True)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
