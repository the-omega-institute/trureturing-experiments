"""Identify the chosen graph source among all laws with the same four weights."""
from collections import Counter
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from runpy import run_path

model = run_path(str(Path(__file__).with_name("ch39-relative-chord-word-source.py")))
local = model["shell"]["local"]


def assignments(hist):
    # Only h=342 has multiplicity one modulo three. The 234 pairs with h=1
    # fit only branch 1. The 168 pairs with h=0 then fit only branch 0.
    assert hist[342] == 1 and all(n % 3 == 0 for h, n in hist.items() if h != 342)
    assert hist[1] == 234 and hist[0] == 168
    free = set(hist)-{342, 1, 0}
    for fourth in free:
        if hist[fourth] != 6:
            continue
        others = sorted(free-{fourth})
        for k in range(len(others)+1):
            for extra in combinations(others, k):
                if sum(hist[h] for h in extra) == 312-234:
                    law = dict.fromkeys(others, 2)
                    law.update(dict.fromkeys(extra, 1))
                    law.update({342:0, 0:0, 1:1, fourth:3})
                    yield law


def verify():
    sites, _ = local["lifted_sites"]()
    lap, _ = local["lifted_laplacian"]()
    second = lap*lap
    hist = Counter(int(h) for row in second.tolist() for h in row)
    labels = sorted(hist)
    target = {h:(0 if h in (0,342) else 1 if 0 < h <= 6 else 3 if h == -3 else 2)
              for h in labels}
    laws = list(assignments(hist))
    assert len(laws) == 1110 and laws.count(target) == 1
    words = model["chord_words"]
    port = [v[0]**2 for v in words]
    unit = sum(v*v for v in words[0])
    assert unit == 5**17 and len(set(port)) == 4

    rivals = [law for law in laws if law != target]
    # A conditional reading at a specified h identifies its branch exactly.
    changed = {h:sum(1 << i for i, law in enumerate(rivals) if law[h] != target[h]) for h in labels}
    need = (1 << len(rivals))-1
    def identifies(picks):
        covered = 0
        for h in picks:
            covered |= changed[h]
        return covered == need
    assert not any(identifies(picks) for picks in combinations(labels, 4))
    assert [picks for picks in combinations(labels, 5) if identifies(picks)] == [(-3, 2, 4, 5, 6)]

    # Three equally frequent values use three different branches. With any
    # fixed binary tag at least two of them coincide and their laws can swap.
    twins = (-21, -3, 5)
    assert {hist[h] for h in twins} == {6} and len({target[h] for h in twins}) == 3
    tag = dict.fromkeys(labels, 0)
    tag.update(dict.fromkeys((2,4,5,6), 2))
    tag[-3] = 1
    def score(law, observable):
        return sum(n*observable[h]*port[law[h]] for h,n in hist.items())
    denominator = 2*625*unit
    graph = score(target, tag)
    gap = min(abs(score(law, tag)-graph) for law in rivals)
    assert graph == 44987627238690 and gap == 822456960
    first = {h:h for h in labels}
    second_moment = {h:h*h for h in labels}
    assert sum(score(law, first) == score(target, first) for law in laws) == 2
    assert sum(score(law, second_moment) == score(target, second_moment) for law in laws) == 1
    distance = local["geometry"]["squared"]
    chords = 0
    for i,a in enumerate(sites):
        for j,b in enumerate(sites):
            branch = sum(distance(a,b) > cut for cut in (152,632,1376))
            chords += tag[int(second[i,j])]*port[branch]
    assert chords == 45006503360706
    print("1110 compatible assignments; minimum five conditional values or a three-valued tag")
    print("normalized graph score", Fraction(graph, denominator),
          "minimum separation", Fraction(gap, denominator),
          "geometric excess", Fraction(chords-graph, denominator))


if __name__ == "__main__":
    verify()
