"""Exact charge-speed and Fourier coefficients of the original shell's rotor lift."""
from collections import defaultdict, deque
from fractions import Fraction
from pathlib import Path
from runpy import run_path

local = run_path(str(Path(__file__).with_name("ch39-three-step-zero-readout.py")))


def charged_shell():
    nodes, edges = local["shell"]["geometry"]()
    adj = [[] for _ in nodes]
    for i, j, amplitude, charge in edges:
        adj[i].append((j, amplitude, charge))
        adj[j].append((i, amplitude, -charge))
    return nodes, adj, nodes.index((33,33,75))


def returns(adj, root):
    result = []
    paths = {(root,0): Fraction(1)}
    for _ in range(4):
        extended = defaultdict(Fraction)
        for (i, winding), amplitude in paths.items():
            for j, weight, charge in adj[i]:
                extended[j,winding+charge] += amplitude*weight
        paths = extended
        result.append({wind: value for (j,wind),value in paths.items() if j == root})
    return result


def reachable(edges, root):
    seen = {root}
    while True:
        extended = seen | {j for i,j,_ in edges if i in seen}
        if extended == seen:
            return seen
        seen = extended


def verify():
    nodes, adj, root = charged_shell()
    edges = [(i,j,charge) for i,row in enumerate(adj) for j,_,charge in row]
    assert len(edges) == 576
    potential = [0]*len(nodes)
    for _ in nodes:
        for i,j,charge in edges:
            potential[j] = max(potential[j], potential[i]+2*charge-1)
    assert all(2*charge-1 <= potential[j]-potential[i] for i,j,charge in edges)
    tight = [(i,j,q) for i,j,q in edges if 2*q-1 == potential[j]-potential[i]]
    core = reachable(tight,root) & reachable([(j,i,q) for i,j,q in tight],root)
    parity, queue = {root:0}, deque([root])
    while queue:
        i = queue.popleft()
        for u,j,q in tight:
            if u != i or j not in core:
                continue
            required = (parity[i]+q)%2
            if j not in parity:
                parity[j] = required
                queue.append(j)
            assert parity[j] == required
    assert len(core) == len(parity) == 18
    r = returns(adj,root)
    assert not r[0] and set(r[1]) == {0}
    for row in r:
        assert all(value > 0 and row[-wind] == value for wind,value in row.items())
    assert set(r[2]) == {-1,0,1} and set(r[3]) == {-2,-1,0,1,2}
    assert r[2][1] == Fraction(3517425,84684693323776)
    assert r[3][2] == Fraction(531441,4736338482895050304)
    print("Closed-walk minimum charged hops 2|k|+(|k| mod 2), tight-component size", len(core))
    print("Three-step and four-step extremal winding coefficients",r[2][1],r[3][2])


if __name__ == "__main__":
    verify()
