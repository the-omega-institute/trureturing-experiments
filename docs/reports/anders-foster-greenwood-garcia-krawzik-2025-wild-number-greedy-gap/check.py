"""Wild number and every legal run of Algorithm 5.10 (Anders-Foster-Greenwood-Garcia-Krawzik, arXiv:2508.06711v1).
Graph: Example 5.12 wheel, and vertex amalgamations of r copies glued in a chain."""
import itertools, sys
from functools import lru_cache

def comps(nv, edges):
    p = list(range(nv))
    def f(x):
        while p[x] != x:
            p[x] = p[p[x]]; x = p[x]
        return x
    for u, v in edges:
        a, b = f(u), f(v)
        if a != b: p[a] = b
    return len({f(x) for x in range(nv)})

def wheel():
    E = {1: [(0,1),(0,3),(0,4)], 2: [(0,2),(0,6),(4,5),(6,7)], 3: [(0,5),(0,7),(1,2),(1,7),(5,6)], 4: [(2,3),(3,4)]}
    return 8, [(u, v, c) for c, es in E.items() for (u, v) in es]

def chain(r):
    nv, es = wheel()
    out = []; off = 0; N = 0
    for k in range(r):
        # copy k uses vertices off..off+7; its vertex 0 is identified with the previous copy's vertex 4
        m = {}
        for x in range(8):
            m[x] = off + x
        if k > 0:
            m[0] = prev4
        for (u, v, c) in es:
            out.append((m[u], m[v], c))
        prev4 = m[4]
        off += 8
    used = sorted({x for (u, v, c) in out for x in (u, v)})
    re = {x: i for i, x in enumerate(used)}
    return len(used), [(re[u], re[v], c) for (u, v, c) in out]

def kappa(nv, es, W, colors):
    return {c: comps(nv, [(u, v) for (u, v, cc) in es if cc == c] + [(es[i][0], es[i][1]) for i in W]) for c in colors}

def color_connected(nv, es, W, colors):
    return all(k == 1 for k in kappa(nv, es, W, colors).values())

def wild_number(nv, es, colors):
    for k in range(len(es) + 1):
        for W in itertools.combinations(range(len(es)), k):
            if color_connected(nv, es, W, colors):
                return k

def dip(nv, es, W, e, colors, base):
    after = kappa(nv, es, tuple(W) + (e,), colors)
    return sum(base[c] - after[c] for c in colors)

def all_runs(nv, es, colors):
    """Every output size of Algorithm 5.10 over all tie-breaking choices. Edges of G/W are the original edges
    (a wild edge contracts its endpoints; an edge inside one class has dip number 0)."""
    sizes = set()
    seen = {}
    def rec(W):
        key = frozenset(W)
        if key in seen: return
        seen[key] = True
        if color_connected(nv, es, W, colors):
            sizes.add(len(W)); return
        base = kappa(nv, es, W, colors)
        cand = [e for e in range(len(es)) if e not in W]
        d = {e: dip(nv, es, W, e, colors, base) for e in cand}
        top = max(d.values())
        best = [e for e in cand if d[e] == top]
        def potential(e):
            W2 = tuple(W) + (e,)
            b2 = kappa(nv, es, W2, colors)
            return tuple(sorted((dip(nv, es, W2, f, colors, b2) if f not in W2 else 0 for f in range(len(es))), reverse=True))
        pot = {e: potential(e) for e in best}
        pm = max(pot.values())
        for e in best:
            if pot[e] == pm:
                rec(tuple(W) + (e,))
    rec(())
    return sizes

colors = [1, 2, 3, 4]
nv, es = wheel()
print("wheel: n=%d m=%d wld=%d greedy_sizes=%s" % (nv, len(es), wild_number(nv, es, colors), sorted(all_runs(nv, es, colors))))
for r in (2, 3):
    nv, es = chain(r)
    print("chain r=%d: n=%d m=%d greedy_sizes=%s" % (r, nv, len(es), sorted(all_runs(nv, es, colors))))

# Fidelity: reproduce Example 5.11 (Table 1 dips, Table 2 potentials, wld = 3) and check that the Example 5.12 order is a legal run.
def ex511():
    idx = {f"v{i}": i - 1 for i in range(1, 7)}
    E = {1: ["v1v3", "v1v4", "v3v6"], 2: ["v2v3", "v4v6", "v5v6"], 3: ["v1v2", "v2v5", "v4v5"]}
    es = []
    for c, l in E.items():
        for s in l:
            a, b = s[:2], s[2:]
            es.append((idx[a], idx[b], c))
    return 6, es
nv, es = ex511(); cols = [1, 2, 3]
name = lambda e: "v%dv%d" % (es[e][0] + 1, es[e][1] + 1)
base = kappa(nv, es, (), cols)
print("Ex5.11 dips:", {name(e): dip(nv, es, (), e, cols, base) for e in range(len(es))})
def potential(nv, es, W, e, cols):
    W2 = tuple(W) + (e,); b2 = kappa(nv, es, W2, cols)
    return tuple(sorted((dip(nv, es, W2, f, cols, b2) if f not in W2 else 0 for f in range(len(es))), reverse=True))
for e in range(len(es)):
    if dip(nv, es, (), e, cols, base) == 2:
        print("  potential", name(e), potential(nv, es, (), e, cols))
print("Ex5.11 wld:", wild_number(nv, es, cols), "greedy sizes:", sorted(all_runs(nv, es, cols)))

def legal_step(nv, es, W, e, cols):
    base = kappa(nv, es, W, cols)
    cand = [f for f in range(len(es)) if f not in W]
    d = {f: dip(nv, es, W, f, cols, base) for f in cand}
    top = max(d.values())
    if d.get(e) != top: return False
    pm = max(potential(nv, es, W, f, cols) for f in cand if d[f] == top)
    return potential(nv, es, W, e, cols) == pm
nv, es = wheel(); cols = [1, 2, 3, 4]
order = [(1, 7), (2, 3), (4, 5), (5, 6), (6, 7), (0, 4)]
pos = {(min(u, v), max(u, v)): i for i, (u, v, c) in enumerate(es)}
W = ()
ok = True
for (u, v) in order:
    e = pos[(min(u, v), max(u, v))]
    ok = ok and legal_step(nv, es, W, e, cols); W = W + (e,)
print("Ex5.12 printed order is a legal run:", ok, "and color-connected at the end:", color_connected(nv, es, W, cols))
