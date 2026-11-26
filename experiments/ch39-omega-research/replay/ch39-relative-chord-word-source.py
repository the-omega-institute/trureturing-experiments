"""Find the least number of short port laws selected by the mutual chord."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import arb, fmpz_mat

shell = run_path(str(Path(__file__).with_name("ch39-shell-selected-short-words.py")))
short, events, q = (shell[k] for k in ("short", "events", "q"))


chord_words = ((-527319, 246668, 336224, 557658), (-526967, 246408, 335994, 558244),
               (527281, 246672, 336138, 557744), (-527151, 246378, 336544, 557752))


def chord_histogram():
    local = shell["local"]
    sites, _ = local["lifted_sites"]()
    squared = local["geometry"]["squared"]
    histogram = Counter(squared(a, b) for a, b in product(sites, repeat=2))
    assert len(histogram) == 31 and sum(histogram.values()) == 625
    assert histogram[0] == 43 and all(n % 6 == 0 for d, n in histogram.items() if d)
    return histogram


def admits_two(bounds, i, j):
    # Intersect four closed half-lines with the whole real mixing segment.
    low, high = Fraction(0), Fraction(1)
    for test in bounds:
        for edge, sign in ((0, 1), (1, -1)):
            a, b = sign*int(test[i, edge]), sign*int(test[j, edge])
            if a > b:
                high = min(high, Fraction(-b, a-b))
            elif a < b:
                low = max(low, Fraction(-b, a-b))
            elif b > 0:
                return False
    return low <= high


def three_word_candidates(bounds, histogram):
    # The zero chord is indivisible. Relax every other multiplicity into
    # independent blocks of six, giving a superset of all possible weights.
    blocks = (625-histogram[0])//6
    triples = []
    for a in range(blocks+1):
        for b in range(blocks-a+1):
            for special in range(3):
                counts = [6*a, 6*b, 6*(blocks-a-b)]
                counts[special] += histogram[0]
                if min(counts) > 0:
                    triples.append(counts)
    weights = np.array(triples, dtype=np.int64)
    assert int(np.max(abs(bounds)))*625 < 2**63
    for ids in combinations(range(bounds.shape[1]), 3):
        tests = bounds[:, ids, :].transpose(1, 0, 2).reshape(3, 4)
        if np.any(tests[:, ::2].min(axis=0) > 0) or np.any(tests[:, 1::2].max(axis=0) < 0):
            continue
        sums = weights @ tests
        keep = np.all(sums[:, ::2] <= 0, axis=1) & np.all(sums[:, 1::2] >= 0, axis=1)
        for counts in weights[keep]:
            yield ids, tuple(map(int, counts))


def weighted_ports(ports, counts):
    return [sum(n*x for n, x in zip(counts, column, strict=True)) for column in zip(*ports, strict=True)]


def batch_masses(port):
    weights = events["orbit_weights_from_ports"](33, port)
    masses, _ = events["event_masses"](weights)
    masses["ref"] -= weights[18, 6, 3]
    return masses


def batch_readout(tables, counts):
    def mass(key):
        return sum(n*table[key] for n, table in zip(counts, tables, strict=True))
    reference = mass("ref")
    return [Fraction(mass(name), reference) for name in events["names"]]


def source_tests(mass, tables, counts, ref, ref_y, logs, constraints):
    dx = [arb(mass[0])/v-q(ref[0]/p) for v, p in zip(mass[1:], ref[1:], strict=True)]
    tests = [(center+sum(a*x for a, x in zip(row[:3], dx, strict=True)), row[3:], margin)
             for row, center, margin in constraints]

    def violates(ys):
        dy = [y-r for y, r in zip(ys, ref_y, strict=True)]
        return any(abs(offset+sum(a*x for a, x in zip(row, dy, strict=True))) > margin
                   for offset, row, margin in tests)

    if violates(logs):
        return ("whole box", "whole box")
    laws = (short["table_readout"](mass), batch_readout(tables, counts))
    result = []
    for ys in laws:
        if any(not events["low"][k] <= v <= events["high"][k]
               for k, v in zip(events["names"], ys, strict=True)):
            result.append("Yukawa")
        elif violates([q(v).log() for v in ys]):
            result.append("gauge")
        else:
            result.append("unresolved")
    return tuple(result)


def exclude_three(ports, histogram):
    ref, ref_y, logs, constraints, _, _ = shell["scale_constraints"](q("26.81672869387187"))
    bounds = shell["affine_port_bounds"](ports, ref, ref_y, logs, constraints)
    assert not any(admits_two(bounds, i, j) for i, j in combinations(range(len(ports)), 2))
    tables = [batch_masses(row) for row in ports]
    excluded = Counter()
    for ids, counts in three_word_candidates(bounds, histogram):
        mass = weighted_ports([ports[i] for i in ids], counts)
        outcome = source_tests(mass, [tables[i] for i in ids], counts, ref, ref_y, logs, constraints)
        if outcome == ("whole box", "whole box"):
            excluded["whole box"] += 1
        else:
            excluded.update(zip(("redraw", "batch"), outcome, strict=True))
    assert excluded == {"whole box": 1438, ("redraw", "Yukawa"): 26, ("redraw", "gauge"): 73,
                        ("batch", "Yukawa"): 19, ("batch", "gauge"): 80}
    print("No continuous two-word mixture; all chord assignments to three words excluded.")


def pointer_gap(cuts, histogram):
    local = shell["local"]
    sites, _ = local["lifted_sites"]()
    squared = local["geometry"]["squared"]
    def band(a, b):
        return sum(squared(a, b) > d for d in cuts)
    increments = {band(c, b)-band(a, b) for a, b, c in product(sites, repeat=3)
                  if 0 < squared(a, c) <= 152}
    assert increments == {-1, 0, 1}
    lap, _ = local["lifted_laplacian"]()
    size = len(sites)
    # On constants this matrix equals size*I, and on their complement
    # it equals size*(2*L-I). Sylvester's criterion gives the spectral gap.
    matrix = [[size*(2*lap[i,j]-int(i == j))+2 for j in range(size)] for i in range(size)]
    assert all(fmpz_mat([r[:k] for r in matrix[:k]]).det() > 0 for k in range(1, size+1))
    pinned_fraction = Fraction(histogram[0], size**2)
    gap = pinned_fraction/3
    assert gap == Fraction(43, 1875)
    print("unit-edge four-state pointer: spectral gap greater than", gap)


def verify():
    histogram = chord_histogram()
    words = shell["words"]
    assert all(not words["hits"](d) for d in range(16))
    library = [(v, d) for d in (16, 17) for v, _ in words["hits"](d)]
    ports = [[x*x*5**(17-d) for x in v] for v, d in library]
    assert len(ports) == len(set(map(tuple, ports))) == 43
    assert all(sum(row) == 5**17 for row in ports)
    exclude_three(ports, histogram)

    cuts = (152, 632, 1376)
    pointer_gap(cuts, histogram)
    counts = [sum(n for d, n in histogram.items() if sum(d > u for u in cuts) == k)
              for k in range(4)]
    assert counts == [169, 312, 138, 6]
    rows = [ports[library.index((v, 17))] for v in chord_words]
    mass = weighted_ports(rows, counts)
    state = short["common_state"](mass)
    batch = batch_readout([batch_masses(v) for v in rows], counts)
    assert all(events["low"][k] < v < events["high"][k] for k, v in zip(events["names"], batch, strict=True))
    scale = q("26.81666203239079")
    for name, initial, bound in (("redraw", state, "2.775e-7"),
                                 ("batch", state[:3]+[q(v) for v in batch], "2.851e-7")):
        errors = short["residuals"](initial, scale)
        assert max(errors) < q(bound) < q("3.18e-7")
        print(name, "four chord bands", counts, "gauge bounds", *(v.upper().str(10) for v in errors))


if __name__ == "__main__":
    verify()
