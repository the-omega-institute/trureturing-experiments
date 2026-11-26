"""十九点壳的 Dirac 物质最小维数：按谱配对给出与裸质量无关的阈值差上界。"""
from collections import Counter
from fractions import Fraction
from itertools import product
from math import prod
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

crossing = run_path(str(Path(__file__).with_name("ch39-unit-seventy-five-crossing.py")))


def dimension(label):
    return prod(Fraction(j-i+sum(label[i:j]), j-i)
                for i in range(4) for j in range(i+1, 5))


def projected_indices(label):
    # 60 times -(T_Y - 3 T_W + 2 T_C) per charge block.
    cells = [(i, j) for i in range(4) for j in range(sum(label[i:]))]
    entries, counts, weights = {}, [0]*5, Counter()
    def fill(k):
        if k == len(cells):
            q = -2*sum(counts[:3])+3*sum(counts[3:])
            weights[abs(q)] += (-q*q + 45*(counts[3]-counts[4])**2
                                - 30*(counts[0]-counts[1])**2)
            return 1
        i, j = cells[k]
        lower = max(entries.get((i, j-1), 0), entries.get((i-1, j), -1)+1)
        total = 0
        for value in range(lower, 5):
            entries[i, j] = value
            counts[value] += 1
            total += fill(k+1)
            counts[value] -= 1
        return total
    assert fill(0) == dimension(label)
    assert sum(weights.values()) == 0
    return {q: w for q, w in weights.items() if w}


def bound(weights, spectra):
    # Pair the positive and negative weights; any ordering yields a valid bound.
    positive, negative = [], []
    for q, w in weights.items():
        side = positive if w > 0 else negative
        side.extend((m, abs(w)) for m in spectra[q])
    positive.sort(key=lambda pair: float(pair[0]))
    negative.sort(key=lambda pair: float(pair[0]))
    total = arb(0)
    i = j = 0
    used_pos = used_neg = 0
    while i < len(positive):
        a, p = positive[i]
        b, n = negative[j]
        take = min(p-used_pos, n-used_neg)
        total += take * (b/a).log().max(arb(0))
        used_pos += take
        used_neg += take
        if used_pos == p:
            i += 1
            used_pos = 0
        if used_neg == n:
            j += 1
            used_neg = 0
    assert j == len(negative)
    return total/(180*arb.pi())


def verify():
    ctx.prec = 224
    # Weyl dimension increases in each label. Every coordinate >=5 exceeds 75.
    assert all(dimension(tuple(5*int(i == j) for i in range(4))) > 75 for j in range(4))
    labels = sorted((int(dimension(l)), l) for l in product(range(5), repeat=4)
                    if l >= l[::-1] and dimension(l) <= 75)
    assert labels[-1] == (75, (0, 1, 1, 0))
    smaller = [(d, l, projected_indices(l)) for d, l in labels[:-1]]
    assert [d for d, l, w in smaller] == [1, 5, 10, 15, 24, 35, 40, 45, 50, 70, 70]
    assert all(w == projected_indices(l[::-1]) for d, l, w in smaller)
    expected = {q: 60*(-y+3*w-2*c) for q,(y,w,c) in crossing["index_table"]().items()}
    assert projected_indices(labels[-1][1]) == expected
    charges = set().union(*(w for d, l, w in smaller))
    assert 0 not in charges
    spectra = crossing["shell_spectra"](charges)
    assert all(len(row) == 19 for row in spectra.values())
    costs = []
    for d, l, w in smaller:
        limit = bound(w, spectra)
        print(d, l, limit.upper().str(10))
        costs.append((d, int((limit.upper()*10**6).ceil().unique_fmpz())))
    # Integer upper bounds in millionths. Each direct summand may be repeated.
    best = [0]
    for budget in range(1, 76):
        best.append(max([best[-1]] + [best[budget-d]+cost for d, cost in costs if d <= budget]))
    assert best[74] <= 2560213 and best[75] <= 2599427
    print("Bounds in millionths for dimensions at most 74, 75", best[74], best[75])
    em, weak = arb("127.930"), arb(".23129")
    initial = (3*(1-weak)*em/5, weak*em, 1/arb(".1180"))
    c_light = initial[0]-3*initial[1]+2*initial[2]
    c_light += (arb("1e17")/arb("91.1876")).log()/(5*arb.pi())
    assert c_light < -arb("10.6083")
    assert c_light + arb(best[75])/10**6 < -arb("8.0088")
    assert c_light + arb(best[74])/10**6 < -arb("8.0480")


if __name__ == "__main__":
    verify()
