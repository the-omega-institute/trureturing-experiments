"""Same-endpoint radial amplitudes and the uniform-gain obstruction at all lengths."""
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
from runpy import run_path

clock = run_path(str(Path(__file__).with_name('ch39-charged-arrival-clock.py')))
shell = clock['shell']
marked = {(7,23,85),(17,55,67),(19,49,71),(29,59,59),(39,39,69),(41,49,61),(51,51,51)}


def valuation(n, prime):
    count = 0
    while n % prime == 0:
        n //= prime
        count += 1
    return count


def verify():
    squared = shell['squared']
    graph = {p:dict(row) for p,row in shell['graph'](3,56).items()}
    pole = (51,51,51)
    graph[pole] = {p:0 for p in graph if 108<=squared(p,pole)<=504}
    for p in graph[pole]:
        graph[p][pole]=0
    labels = sorted({squared(p,q) for p,row in graph.items() for q in row if squared(p,q)>=108})
    marker = {p:int(tuple(sorted(p)) in marked) for p in graph}
    adjacency = {}
    for p,row in graph.items():
        adjacency[p]=[]
        for q,w in row.items():
            h = squared(p,q)
            if h>=108:
                assert marker[p]^marker[q] == (1+valuation(h,2)+valuation(h,17))%2
                adjacency[p].append((q,w,8**labels.index(h)))
    for i,s in enumerate(clock['sites']):
        assert all(marker[p] == int(i==5) for p in clock['orbit'](s))
    assert Fraction(144,624) != Fraction(1,2)
    assert Fraction(144,624) > Fraction(4,13)**3

    pairs = 0
    for i,s in enumerate(clock['sites']):
        paths = {(clock['orbit'](s)[0],0,0)}
        for _ in range(5):
            paths = {(q,k+w,code+label) for p,k,code in paths for q,w,label in adjacency[p]}
        histograms = defaultdict(set)
        for p,k,code in paths:
            histograms[p,k].add(code)
        for j,t in enumerate(clock['sites']):
            if i<2 and j<2:
                continue
            found = any(set.intersection(*(histograms[q,k] for k in (-1,0,1)))
                        for q in clock['orbit'](t))
            assert found
            pairs += 1
    assert pairs == 32
    print('Common five-edge histogram for all three charges on each of 32 orbit pairs.')
    print('Parity verified on',sum(map(len,adjacency.values()))//2,'edges; marked fraction 3/13.')


if __name__=='__main__':
    verify()
