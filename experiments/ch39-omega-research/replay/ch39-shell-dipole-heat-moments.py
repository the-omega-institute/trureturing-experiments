"""十九点无磁单位壳上，端口热矩、有效电阻和单条内部边的可辨识阶数。"""
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
from runpy import run_path

import sympy as sp


def verify():
    geometry = run_path(str(Path(__file__).with_name("ch39-arc-relative-determinant.py")))
    sites = geometry["sites"]
    identity = sp.eye(len(sites))
    lap = sp.zeros(len(sites))
    edges = [(i, j) for i, j, _, _ in geometry["edges"]]
    for i, j in edges:
        edge = identity[:, i] - identity[:, j]
        lap += edge * edge.T
    projector = sp.ones(len(sites)) / len(sites)
    green = (lap + projector).inv()
    assert lap * green == identity - projector
    moment_groups = defaultdict(set)
    first_visible = Counter()
    saturated = 0
    for a, b in combinations(range(len(sites)), 2):
        source = identity[:, a] - identity[:, b]
        powers = [source]
        for _ in range(1, len(sites)):
            powers.append(lap * powers[-1])
        p = source.dot(powers[1]) / 2
        q = powers[1].dot(powers[1]) / 2
        potential = green[:, a] - green[:, b]
        response = source.dot(potential) / 2
        bound = max(1 / p, 2 * p / (4 * p**2 - q))
        assert p**2 <= q < 4 * p**2 and response >= bound
        assert (response == bound) == (q == p**2)
        saturated += int(response == bound)
        moment_groups[p, q].add(response)

        for i, j in edges:
            if i in (a, b) or j in (a, b):
                continue
            order = next((k for k, vector in enumerate(powers[:5]) if vector[i] != vector[j]), None)
            first_visible[order] += 1
            if order is None:
                assert potential[i] == potential[j]
                # 全部高阶矩相同需验证端口的整个 Krylov 空间。
                assert all(vector[i] == vector[j] for vector in powers[5:])
            else:
                assert potential[i] != potential[j]
    assert saturated == 12 and len(moment_groups) == 15
    assert all(len(responses) == 1 for responses in moment_groups.values())
    assert first_visible == {1: 6612, 2: 5940, 3: 1812, 4: 144, None: 3852}
    print("单位壳 171 个端口的两矩类型、达到下界数", len(moment_groups), saturated)
    print("内部边首次可辨的热矩阶数 3、5、7、9、全阶不可辨",
          [first_visible[k] for k in (1, 2, 3, 4, None)])

    a, b, i, j = [sites.index(site) for site in
                  [(23, 17, 7), (25, 11, 11), (11, 11, 25), (23, 7, 17)]]
    assert (i, j) in edges
    source = identity[:, a] - identity[:, b]
    edge = identity[:, i] - identity[:, j]
    changed = lap + edge * edge.T
    before, after = source, source
    for k in range(1, 10):
        before, after = lap * before, changed * after
        assert source.dot(after-before) / 2 == (200 if k == 9 else 0)
    old = (source.T * green * source)[0]
    new = (source.T * (changed + projector).inv() * source)[0]
    assert old - new == sp.Rational(6800, 3384036987167)
    print("指定端口单位电流电压差的减少量，第九热矩的增加量", old-new, 200)


if __name__ == "__main__":
    verify()
