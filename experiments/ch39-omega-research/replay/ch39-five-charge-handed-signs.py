"""五个荷的三阶迹手性差在实际零点上的符号模式和自然密度。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


CHARGES = (1, 2, 3, 4, 6)


def certify():
    root = Path(__file__).resolve().parent
    exact = run_path(str(root / "ch39-mincut-real-triangle-tomography.py"))
    ctx.prec = 160
    units = [exact["unit"](label) for label in ((510, 1), (1173, -2))]
    coefficients = [tuple(-36*(unit**q).imag for unit in units) for q in CHARGES]
    n = 512
    nodes = [2*arb.pi()*i/n for i in range(n+1)]
    cells = [l.union(r) for l, r in zip(nodes, nodes[1:])]
    first, second = [], []
    for q, (a, b) in zip(CHARGES, coefficients, strict=True):
        assert not a.contains(0) and not b.contains(0)
        waves = [(q*t).sin() for t in cells]
        first.append([a*s for s in waves])
        second.append([b*s for s in waves])
    counts = [0]*32
    for i in range(n):
        for j in range(n):
            mask = 0
            for k in range(5):
                value = first[k][i] + second[k][j]
                if value.contains(0):
                    break
                if value > 0:
                    mask |= 1 << k
            else:
                counts[mask] += 1
    assert min(counts)*100 > n*n
    print("二维等面积整格数、单模式最小占格数、最小密度下界",
          n*n, min(counts), (arb(min(counts))/(n*n)).str(15))
    # 始终按实际编号递增。球包含零时不作符号判定。
    geo = exact["graph"]
    first_seen = {}
    for number, zero in enumerate(acb.zeta_zeros(1, 250), 1):
        phases = geo["phases"](zero.imag)
        pattern = 0
        for k, (q, (a, b)) in enumerate(zip(CHARGES, coefficients, strict=True)):
            value = a*(phases[724]**q).imag + b*(phases[804]**q).imag
            assert not value.contains(0)
            pattern |= int(value > 0) << k
        first_seen.setdefault(pattern, number)
    assert len(first_seen) == 32, set(range(32))-first_seen.keys()
    print("按五位掩码 0..31 排列的首枚实际编号", [first_seen[i] for i in range(32)])
    return first_seen


if __name__ == "__main__":
    certify()
