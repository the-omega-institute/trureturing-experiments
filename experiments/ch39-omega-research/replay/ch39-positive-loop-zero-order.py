"""Sharp length at which a positive loop score can reverse a high-zero ordering."""
from pathlib import Path
from runpy import run_path

from flint import acb, ctx

shell = run_path(str(Path(__file__).with_name("ch39-refined-mass-obstruction.py")))


def winding_certificate(nodes, edges):
    directed = {(i,j): k for i,j,_,k in edges}
    directed.update({(j,i): -k for i,j,_,k in edges})
    potential = [0]*len(nodes)
    for _ in nodes:
        previous = potential.copy()
        for (i,j), k in directed.items():
            potential[j] = max(potential[j], potential[i]+2*k-1)
        if previous == potential:
            break
    assert all(2*k-1 <= potential[j]-potential[i] for (i,j), k in directed.items())
    tight = [(i,j) for (i,j), k in directed.items() if 2*k-1 == potential[j]-potential[i]]
    reachable = [{i} for i in range(len(nodes))]
    for _ in nodes:
        for i,j in tight:
            reachable[i].update(reachable[j])
    cyclic = [(i,j) for i,j in tight if i in reachable[j]]
    assert len(cyclic) == 24
    color = {}
    for source,_ in cyclic:
        if source in color:
            continue
        color[source] = 0
        pending = [source]
        while pending:
            i = pending.pop()
            for u,j in cyclic:
                if u == i and j not in color:
                    color[j] = (color[i]+1)%4
                    pending.append(j)
    assert all(color[j] == (color[i]+1)%4 for i,j in cyclic)
    square = [(21,51,69), (33,33,75), (51,21,69), (39,39,69), (21,51,69)]
    triangle = [square[0], square[1], square[3], square[0]]
    for path, winding in [(square, 2), (triangle, 1)]:
        idx = [nodes.index(p) for p in path]
        assert sum(directed[i,j] for i,j in zip(idx,idx[1:])) == winding


def verify():
    ctx.prec = 192
    nodes, edges = shell["geometry"]()
    winding_certificate(nodes, edges)
    phases = []
    for n in range(11485, 11492):
        rho = acb.zeta_zero(n)
        phases.append((rho-1)/rho)
    # Ball multiplications enumerate every harmonic needed, with no angle rounding.
    powers = [acb(1) for _ in phases]
    first_drop = [None]*(len(phases)-1)
    for k in range(1, 34996):
        powers = [a*b for a,b in zip(powers, phases, strict=True)]
        for i, (lo, hi) in enumerate(zip(powers, powers[1:])):
            if first_drop[i] is None:
                if hi.real < lo.real:
                    first_drop[i] = k
                else:
                    assert hi.real > lo.real
    assert first_drop == [34983,34986,34988,34991,34993,34995]
    lengths = [2*k + k%2 for k in first_drop]
    assert lengths[0] == 69967
    print("Adjacent pair first reversed winding", first_drop)
    print("Minimum closed walk lengths", lengths)


if __name__ == "__main__":
    verify()
