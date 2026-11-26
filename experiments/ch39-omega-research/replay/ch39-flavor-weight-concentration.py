"""两个固定平均跃迁概率对单枚零点的味概率及离散选权熵的限制。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

joint = run_path(str(Path(__file__).with_name("ch39-three-window-two-probability-cp.py")))


def binary_entropy(p):
    return -(p*p.log()+(1-p)*(1-p).log())/arb(2).log()


def verify():
    centers = [joint["inner_centers"][0], (".1225", ".01581205"), joint["inner_centers"][1]]
    locations = joint["previous"]["open_windows"](centers)
    jets = joint["observables"]()
    low, high = map(arb, joint["outer_interval"])
    a, b = joint["box_values"](jets, low)[:2]
    x, y = joint["box_values"](jets, high)[:2]
    slope = (b-y)/(a-x)
    intercept = b-slope*a
    up, charm, _ = jets(low.union(high))
    assert up[1] < 0 and charm[2]*up[1]-up[2]*charm[1] < 0

    p = arb(".686803")
    pending = [(low, high)]
    while pending:
        left, right = pending.pop()
        up, charm, _ = joint["box_values"](jets, left.union(right))
        other = -p*up/(1-p)
        chord_gap = intercept+slope*up-charm
        if other < x or other > a or p*chord_gap > intercept-2:
            continue
        assert right-left > arb("1e-12")
        mid = (left+right)/2
        pending.extend(((left, mid), (mid, right)))

    rows = [joint["box_values"](jets, c) for c in locations]
    matrix = arb_mat([[1]*3, [r[0] for r in rows], [r[1] for r in rows]])
    weights = matrix.solve(arb_mat([[1], [0], [2]])).entries()
    assert all(w > 0 for w in weights)
    entropy = -sum(w*w.log() for w in weights)/arb(2).log()
    assert weights[1] > arb(".686762") and entropy < arb(".89689")
    assert binary_entropy(p) > arb(".8968294")
    print("任一单点的味概率上界、三枚真实零点的权重及熵",
          p.str(10), *[w.str(12) for w in weights], entropy.str(12))


if __name__ == "__main__":
    verify()
