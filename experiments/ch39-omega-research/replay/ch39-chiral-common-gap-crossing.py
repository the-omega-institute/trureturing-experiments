"""十九点两步质量的共同三代分界及其一圈尖阈值联合检验。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

ROOT = Path(__file__).parent
source = run_path(str(ROOT / "ch39-fixed-planck-magnetic-window.py"))
EXPECTED = [10864, 11062, 11198, 11325, 11461, 11597, 11661, 11733, 11870, 11934, 11998]


def ordered(values):
    # 两端分别排序，对有简并的谱也包住每个顺序统计量。
    lo = sorted(v.lower() for v in values)
    hi = sorted(v.upper() for v in values)
    return [a.union(b) for a, b in zip(lo, hi)]


def exclude_crossing(model, initial, layers, trim_higgs):
    b = model["beta"]
    d = [b[1] - b[2], b[2] - b[0], b[0] - b[1]]
    def project(w):
        return sum(a * z for a, z in zip(d, w)), w[0] - w[1]

    terms = []
    def add(values, weight):
        proj, ew = project(weight)
        terms.extend((x.log(), proj, ew) for x in values)

    for charge, family in zip(model["charges"], model["weyl"]):
        add(layers[charge][3:], [2 * x for x in family])
    add(layers[3][3:] if trim_higgs else layers[3], model["higgs"])
    add(layers[0], [arb(0), arb(-7), -arb(21) / 2])

    offset, diff = project([2 * arb.pi() * a for a in initial])
    assert not offset.contains(0)
    L = (arb("1e13") / model["mz"]).log().union((arb("1e17") / model["mz"]).log())
    # t 在最小阈值以下时投影恒为 offset。所有新质量在 M_Z 以上给出另一端。
    lo = min(x.lower() for x, _, _ in terms)
    hi = min(x.upper() for x, _, _ in terms) + L.upper()
    pending = [(lo, hi.upper())]
    while pending:
        left, right = pending.pop()
        t = left.union(right)
        residual, implied = offset, diff
        for x, proj, ew in terms:
            exposed = (t - x).max(0)
            residual -= proj * exposed
            implied -= ew * exposed
        if not residual.contains(0) or not (implied / (b[0] - b[1])).overlaps(L):
            continue
        assert right - left > arb('1e-8'), "未能排除交点"
        mid = ((left + right) / 2).mid()
        pending.extend(((left, mid), (mid, right)))


def verify(start=10855, count=1238):
    model, _, _, initial, _ = source["mean"]["inputs"]()
    initial[2] = 1 / arb(".117").union(arb(".119"))
    ctx.prec = 128
    read = source["spectral_supplier"](model)
    ctx.prec = 128
    hits, fixed = {}, []
    zeros = acb.zeta_zeros(start, count)
    assert len(zeros) == count
    for n, zero in enumerate(zeros, start):
        rows = read(zero.imag)
        layers = {q: ordered([x for x, charge in rows if charge == q])
                  for q in (0, *model["charges"])}
        left = right = None
        for q in model["charges"]:
            third, fourth = layers[q][2:4]
            left = third if left is None else left.max(third)
            right = fourth if right is None else right.min(fourth)
        if left < right:
            hits[n] = (left, right)
            if left < 9 < right:
                fixed.append(n)
            else:
                assert 9 < left or right < 9
            for trim in (False, True):
                exclude_crossing(model, initial, layers, trim)
        else:
            assert left > right
    assert list(hits) == [n for n in EXPECTED if start <= n < start + count]
    assert fixed == [n for n in EXPECTED[1:] if start <= n < start + count]
    if 11661 in hits:
        left, right = hits[11661]
        assert all(right - left > b - a for n, (a, b) in hits.items() if n != 11661)
        assert arb('7.49571') < left < arb('7.49572')
        assert arb('10.59934') < right < arb('10.59936')
        print("11661 分界两端", left.str(12), right.str(12))
    print("共同分界", list(hits), "分界 9", fixed, "两种标量计数的交点均已排除")


if __name__ == "__main__":
    verify()
