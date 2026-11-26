"""在第 2 和 30 枚实际零点，对一切正截止证明两个热核导数严格有序。"""
from pathlib import Path
from runpy import run_path

from flint import arb


def verify():
    root = Path(__file__).resolve().parent
    fields = run_path(str(root / "ch39-two-zero-cutoff-gravity.py"))["gravity_fields"]
    exp_difference = run_path(str(root / "ch39-two-zero-cutoff-support.py"))["exponential_difference"]
    _, _, (first, second) = fields(root, (2, 30))
    pairs = []
    for (g, a, z), (other_g, other_a, h) in zip(first[1:], second[1:], strict=True):
        assert other_g == g and (other_a-a).contains(0)
        pairs.append((a, -g, z, h))

    # 差一律取第 30 枚减第 2 枚。相同的无质量项已消去。
    lead = min(range(len(pairs)), key=lambda j: float(pairs[j][3]))
    alpha, spin, _, ref = pairs[lead]
    tail = [(a, s, m) for a, s, z, h in pairs for m in (h, z)]
    del tail[2*lead]
    gap = min((m-ref).lower() for _, _, m in tail)
    start = arb("1e-14")
    assert gap > 0
    rem = (-gap/start).exp()
    assert alpha > rem*sum(abs(a) for a, _, _ in tail)
    assert spin > rem*sum(abs(s) for _, s, _ in tail)

    # 大 x 时 y=1/x。乘以 x 后在 y=0 有普通极限，故包围无上限的所有截止。
    end = arb(".001")
    y = arb(0).union(1/end)
    limits = [(z-h)*(-y*z.union(h)).exp() for _, _, z, h in pairs]
    newton = sum(a*d for (a, _, _, _), d in zip(pairs, limits))
    vacuum = sum(s*d for (_, s, _, _), d in zip(pairs, limits))
    assert newton > 0 and vacuum > 0

    def ordered(x):
        diffs = [exp_difference(x, z, h, ref) for _, _, z, h in pairs]
        return (sum(a*d for (a, _, _, _), d in zip(pairs, diffs)) > 0
                and sum(s*d for (_, s, _, _), d in zip(pairs, diffs)) > 0)

    count = 0
    left = start
    while left < end:
        right = min(2*left, end)
        pending = [(left, right)]
        while pending:
            lo, hi = pending.pop()
            if ordered(lo.union(hi)):
                count += 1
            else:
                assert hi-lo > arb("1e-22"), (lo, hi)
                mid = (lo+hi)/2
                pending.extend(((lo, mid), (mid, hi)))
        left = right
    print("2 与 30 的红外质量间隙、中段闭区间数、紫外两差的下界",
          gap.str(12), count, newton.lower().str(12), vacuum.lower().str(12))


if __name__ == "__main__":
    verify()
