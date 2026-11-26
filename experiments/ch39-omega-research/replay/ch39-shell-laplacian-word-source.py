"""Synthesize the four branch weights using only a quadratic graph readout."""
from collections import Counter, deque
from itertools import product
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import arb

chord = run_path(str(Path(__file__).with_name("ch39-relative-chord-word-source.py")))
shell = chord["shell"]
short, q = shell["short"], shell["q"]


def graph_rules():
    local = shell["local"]
    sites, _ = local["lifted_sites"]()
    lap, _ = local["lifted_laplacian"]()
    matrix = np.array(lap.tolist(), dtype=np.int64)
    size = len(sites)
    second = matrix @ matrix
    most = second.max()
    near = second[second < 0].max()
    least_diagonal = second.diagonal().min()
    def branch(h):
        if h == 0 or h == most:
            return 0
        if 0 < h <= least_diagonal:
            return 1
        if h == near:
            return 3
        return 2
    assert (most, near, least_diagonal) == (342, -3, 6)
    assert Counter(map(branch, second.flat)) == {0:169, 1:312, 2:138, 3:6}
    assert np.count_nonzero(matrix == 0) == 474 and np.all(matrix.diagonal() > 0)
    # The transposition is a graph automorphism fixing the first endpoint.
    a, b, c = (sites.index(v) for v in ((5,29,1), (1,5,29), (5,1,29)))
    perm = list(range(size))
    perm[b], perm[c] = c, b
    assert np.array_equal(matrix[np.ix_(perm, perm)], matrix)
    distance = local["geometry"]["squared"]
    assert (distance(sites[a], sites[b]), distance(sites[a], sites[c])) == (1376, 1568)
    # Product-graph distance to a pair of coincident physical coordinates.
    meeting = {(i,j):0 for i,j in product(range(size), repeat=2) if sites[i] == sites[j]}
    queue = deque(meeting)
    neighbors = [list(np.flatnonzero(row < 0)) for row in matrix]
    while queue:
        i, j = queue.popleft()
        adjacent = [(v,j) for v in neighbors[i]] + [(i,v) for v in neighbors[j]]
        for x in adjacent:
            if x not in meeting:
                meeting[x] = meeting[i,j]+1
                queue.append(x)
    counts = Counter(meeting.values())
    assert counts == {0:43, 1:126, 2:288, 3:144, 4:24}
    return tuple(counts[i] for i in range(5))


def coarse_candidates(bounds, counts):
    prefixes = np.array(list(product(range(43), repeat=len(counts)-2)))
    initial = sum(w*bounds[:, prefixes[:,k], :] for k,w in enumerate(counts[:-2]))
    for i, j in product(range(43), repeat=2):
        end = counts[-2]*bounds[:,i,:]+counts[-1]*bounds[:,j,:]
        keep = np.ones(len(prefixes), dtype=bool)
        for k in range(2):
            keep &= (initial[k,:,0] <= -end[k,0]) & (initial[k,:,1] >= -end[k,1])
        for prefix in prefixes[keep]:
            yield (*map(int, prefix), i, j)


def final_hop_exclusion(mass, ys, beta, mix):
    time = q("26.81668508")
    shortest = shell["shortest"]
    state = [arb(v)/mass[0] for v in mass[1:]] + [q(v) for v in ys]
    for _ in range(128):
        state = shortest["reference"]["advance"](state, time/128, 14)
    target = shell["model"]["target"](11431, time)
    err = [1/(r*u)-target for r,u in zip(shortest["ratios"], state)]
    b = beta[0]/beta[1]
    slack = abs(err[0]-b*err[1])-(1+b)*(arb(2)**-21+q("1.3e-16"))-mix[0]-b*mix[1]
    assert slack > q("1e-6")
    return slack.lower()


def verify():
    counts = graph_rules()
    ref, ref_y, logs, constraints, beta, mix = shell["scale_constraints"](q("26.81672869387187"))
    words = shell["words"]
    library = [(v,d) for d in (16,17) for v,_ in words["hits"](d)]
    ports = [[x*x*5**(17-d) for x in v] for v,d in library]
    bounds = shell["affine_port_bounds"](ports, ref, ref_y, logs, constraints)
    assert len(ports) == 43 and 625*int(np.max(abs(bounds))) < 2**63
    # Any affine polynomial in L is constant on all off-diagonal nonedges.
    hull = np.stack((bounds[:,:,0].min(axis=1), bounds[:,:,1].max(axis=1)), axis=1)
    boxes = 474*bounds + 151*hull[:,None,:]
    assert not np.any(np.all((boxes[:,:,0] <= 0) & (boxes[:,:,1] >= 0), axis=0))
    tables = [chord["batch_masses"](row) for row in ports]
    outcomes = Counter()
    for ids in coarse_candidates(bounds, counts):
        rows = [ports[i] for i in ids]
        mass = chord["weighted_ports"](rows, counts)
        selected = [tables[i] for i in ids]
        result = chord["source_tests"](mass, selected, counts, ref, ref_y, logs, constraints)
        outcomes.update(zip(("redraw", "batch"), result, strict=True))
        if "unresolved" in result:
            assert ids == (29,34,7,31,34) and result == ("unresolved","unresolved")
            for ys in (short["table_readout"](mass), chord["batch_readout"](selected, counts)):
                print("last meeting-time assignment, exclusion margin", final_hop_exclusion(mass, ys, beta, mix))
    assert outcomes == {("redraw","whole box"):439, ("batch","whole box"):439,
                        ("redraw","Yukawa"):8, ("batch","Yukawa"):3,
                        ("redraw","gauge"):31, ("batch","gauge"):36,
                        ("redraw","unresolved"):1, ("batch","unresolved"):1}
    four = (169, 312, 138, 6)
    outcomes = Counter()
    for ids in coarse_candidates(bounds, four):
        mass = chord["weighted_ports"]([ports[i] for i in ids], four)
        result = chord["source_tests"](mass, [tables[i] for i in ids], four, ref, ref_y, logs, constraints)
        outcomes[result] += 1
        if "unresolved" in result:
            assert result == ("unresolved","unresolved")
            assert [library[i] for i in ids] == [(v, 17) for v in chord["chord_words"]]
    assert outcomes == {("whole box","whole box"):76, ("gauge","gauge"):4,
                        ("Yukawa","gauge"):3, ("unresolved","unresolved"):1}
    print("unique four-class assignment in both protocols")
    print("all 43^5 meeting-distance assignments excluded; minimum Laplacian polynomial degree 2")


if __name__ == "__main__":
    verify()
