"""认证第七壳真空核质量简并度的指数和在全非负半轴为负。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify():
    root = Path(__file__).resolve().parent
    ctx.prec = 192
    source = run_path(str(root/"ch39-three-index-free-atoms.py"))
    field = source["broken_fields"](root, (2306888100092,))[2][0]
    assert all(z >= 0 for _, _, z in field)
    grid = [arb(0)] + [arb(2)**(arb(k)/8) for k in range(257)]
    for lo, hi in zip(grid[:-1], grid[1:], strict=True):
        upper = sum(g*(-z*(lo if g>0 else hi)).exp() for g, _, z in field)
        assert upper < arb("-30.47")
    tail = sum(g*(-z*grid[-1]).exp() if g>0 else g for g, _, z in field if g>0 or z==0)
    assert tail < arb("-30.47")
    print("带号质量简并度的指数和在全半轴严格小于 -30.47")


if __name__ == "__main__":
    verify()
