#!/usr/bin/env python3
"""Exact A185766 checks. Standard library only; --b-file allows offline source data.
Coefficients are supplied by OEIS, never fitted to computed terms.
"""
import argparse
import hashlib
from itertools import product
from pathlib import Path
from urllib.request import Request, urlopen

C = (1, 54, 60, -400, -155, 924, -1220, -76, 1725, -642, -412, -88, 100, 40)
PREFIX = (1736, 10592, 75192, 562800, 4430912, 35108016, 278758936,
          2212654240, 17567479016, 139478572544, 1107370142640, 8792101180232,
          69804143650776, 554214418338440, 4400166223484904,
          34935244761367808, 277367986882503912)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--b-file', type=Path)
    args = parser.parse_args()
    if args.b_file:
        raw = args.b_file.read_bytes()
    else:
        req = Request('https://oeis.org/A185766/b185766.txt',
                      headers={'User-Agent': 'A185766 exact checker'})
        with urlopen(req, timeout=60) as response:
            raw = response.read()
    source = dict(tuple(map(int, s.split())) for s in raw.decode().splitlines()
                  if s.strip() and not s.lstrip().startswith('#'))
    assert set(source) == set(range(1, 201))
    print('OEIS b-file sha256:', hashlib.sha256(raw).hexdigest())
    rows = list(product((0, 1), repeat=7))
    def tau(r, s):
        return tuple(r[j] + s[j+1] for j in range(6))
    traces = {(r, s): tau(r, s) for r in rows for s in rows}
    states = [p for p, z in traces.items() if all(z[j] != z[j+1] for j in range(5))]
    idx = {p: i for i, p in enumerate(states)}
    edges = [[idx[s, t] for t in rows if (s, t) in idx
              and all(x != y for x, y in zip(traces[r, s], traces[s, t]))]
             for r, s in states]
    # Independent bit representation, direct 2-row and 3-row enumeration.
    bt = [[tuple(((r >> j) & 1) + ((s >> (j+1)) & 1) for j in range(6))
           for s in range(128)] for r in range(128)]
    ok = [[all(z[j] != z[j+1] for j in range(5)) for z in line] for line in bt]
    bit_states = {(r, s) for r in range(128) for s in range(128) if ok[r][s]}
    bit_edges = set()
    direct_three = 0
    for r in range(128):
        for s in range(128):
            if ok[r][s]:
                for t in range(128):
                    if ok[s][t] and all(bt[r][s][j] != bt[s][t][j] for j in range(6)):
                        direct_three += 1
                        bit_edges.add(((r, s), (s, t)))
    def encode(r):
        return sum(x << j for j, x in enumerate(r))
    encoded = [(encode(r), encode(s)) for r, s in states]
    assert set(encoded) == bit_states
    assert {(encoded[i], encoded[j]) for i, es in enumerate(edges) for j in es} == bit_edges
    assert len(states) == 1736 and len(bit_edges) == direct_three == 10592
    print('States:', len(states), 'Edges:', len(bit_edges))
    print('Independent bit/tuple states and edges: identical')
    print('Direct two-row and three-row counts:', len(bit_states), direct_three)
    examples = [((0,0,0,1,1,0,0),(1,0,1,0,0,1,0),(1,1,0,0,0,1,1)),
                ((0,1,0,1,1,0,0),(1,0,1,0,0,1,0),(0,0,0,1,0,1,1)),
                ((0,1,1,1,0,0,0),(1,0,1,0,1,1,0),(1,0,0,1,1,1,0))]
    for a in examples:
        z = [tau(a[i], a[i+1]) for i in range(2)]
        assert all(z[i][j] != z[i][j+1] for i in range(2) for j in range(5))
        assert all(z[0][j] != z[1][j] for j in range(6))
    print('All 3 displayed source examples satisfy literal trace constraints')
    def step(v):
        return [sum(v[j] for j in es) for es in edges]
    v = [[1] * len(states)]
    for _ in range(199):
        v.append(step(v[-1]))
    counts = {n: sum(v[n-1]) for n in range(1, 201)}
    assert counts == source
    assert tuple(counts[n] for n in range(1, 18)) == PREFIX
    print('All 200 b-file values match exact transfer counts')
    print('First 6 counts:', [counts[n] for n in range(1, 7)])
    horner = [1] * len(states)
    for c in C:
        horner = [x - c for x in step(horner)]
    metrics = []
    for k in range(7):
        residual = [v[k+14][i] - sum(c*v[k+13-j][i] for j, c in enumerate(C))
                    for i in range(len(states))]
        assert residual == horner
        metrics.append((sum(x != 0 for x in residual), sum(residual), max(map(abs, residual))))
        print('Residual k=%d: nonzero=%d sum=%d max_abs=%d' % ((k,) + metrics[-1]))
        horner = step(horner)
    assert metrics[4:] == [(1652, 0, 229168), (984, 0, 5760), (0, 0, 0)]
    assert metrics[3][1] == -7776
    for n in range(19, 201):
        assert counts[n] == sum(c*counts[n-j] for j, c in enumerate(C, 1))
    assert counts[18] - sum(c*counts[18-j] for j, c in enumerate(C, 1)) == -7776
    print('Printed recurrence checked at n=19..200; residual at n=18: -7776')
    print('T^6 q(T)1=0 coordinatewise; scalar sums at k=4,5 vanish')
    print('Direct and Horner evaluations agree at k=0..6')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
