"""检验右弱壳的标记与广义 CP，枚举定向边稳定子。"""
from itertools import combinations, permutations
from pathlib import Path
from runpy import run_path

flavour = run_path(str(Path(__file__).with_name("ch39-shell-plaquette-flavour.py")))
shell = flavour["shell"]


def automorphisms(neighbours):
    n = len(neighbours)
    adjacent = [set(row) for row in neighbours]
    signatures = [(len(row),sorted(len(adjacent[j]) for j in row)) for row in adjacent]
    candidates = [[j for j in range(n) if signatures[j]==sig] for sig in signatures]
    order = sorted(range(n), key=lambda i:(len(candidates[i]),-len(adjacent[i])))
    mapping, used = {}, set()

    def visit(depth):
        if depth == n:
            yield tuple(mapping[i] for i in range(n))
            return
        i = order[depth]
        for j in candidates[i]:
            if j in used or any((k in adjacent[i])!=(v in adjacent[j]) for k,v in mapping.items()):
                continue
            mapping[i] = j
            used.add(j)
            yield from visit(depth+1)
            used.remove(j)
            del mapping[i]

    return list(visit(0))


def verify():
    _,_,_,right = shell["neutral"]["graph_constants"](port_builder=shell["capture"])
    lap,edges = right
    neighbours, _ = flavour["cycles"](lap,edges)
    hubs = [(i,j) for i,j in combinations(range(len(neighbours)),2)
            if len(set(neighbours[i])&set(neighbours[j]))==3]
    assert len(hubs)==1
    a,b = hubs[0]
    middle = sorted(set(neighbours[a])&set(neighbours[b]))
    group = automorphisms(neighbours)
    reverse = []
    actions = set()
    for perm in group:
        labels = tuple(middle.index(perm[i]) for i in middle)
        assert {perm[a], perm[b]}=={a,b}
        swap = int(perm[a]==b)
        parity = sum(labels[i]>labels[j] for i,j in combinations(range(3),2))%2
        actions.add((labels,swap))
        if parity^swap:
            reverse.append(perm)
    assert actions == {(labels,swap) for labels in permutations(range(3)) for swap in (0,1)}

    def protected(marked):
        return any(all(perm[v]==v for edge in marked for v in edge) for perm in reverse)

    assert all(protected([edge]) for edge in edges)
    broken = [pair for pair in combinations(edges,2) if not protected(pair)]
    assert (len(group),len(reverse),len(edges),len(broken))==(96,48,15,48)
    example = tuple((max(a,j),min(a,j)) for j in middle[:2])
    assert not protected(example)
    print("graph automorphisms",len(group),"conjugating",len(reverse),"single marks",len(edges),
          "two-mark choices breaking conjugation",len(broken),"of",len(list(combinations(edges,2))),flush=True)
    print("same-hub edges breaking conjugation",example,flush=True)


if __name__=="__main__":
    verify()
