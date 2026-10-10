# Check of the Cha–Lee conjecture (arXiv:2501.11816v3, Sec. 6.x): for (k,m)-balanced distributions of the
# QFT circuit (n = km qubits, final SWAPs omitted), the canonical allocation pi*(q_i) = p_{ceil(i/m)} is optimal for MS-GC.
import itertools, random, sys
random.seed(1)
def times(n):
    # conventional order: H_0, CP(0,1..n-1), H_1, CP(1,2..), ... ; times start at 1 (0 is the initial migration slot)
    tH = {}; tG = {}; t = 1
    for i in range(n):
        tH[i] = t; t += 1
        for j in range(i+1, n):
            tG[(i, j)] = t; t += 1
    return tH, tG
def latest(T, tstar):
    return max(x for x in ({0} | T) if x <= tstar)
def candidates(n, k, pi, tH):
    return {(q, p, t) for q in range(n) for p in range(k) if p != pi[q] for t in (0, tH[q])}
def feasible_literal(n, k, pi, S, tH, tG):
    for (i, j), ts in tG.items():
        if pi[i] == pi[j]: continue
        si = latest({tH[i]}, ts); sj = latest({tH[j]}, ts)
        if (i, pi[j], si) in S or (j, pi[i], sj) in S: continue          # home coverage, Definition 3
        if any((i, p, si) in S and (j, p, sj) in S for p in range(k) if p not in (pi[i], pi[j])): continue  # joint coverage, Definition 4
        return False
    return True
def setpairs(n, pi, S, tH):
    L = {q: {pi[q]} for q in range(n)}; R = {q: {pi[q]} for q in range(n)}
    for (q, p, t) in S:
        (L if t == 0 else R)[q].add(p)
    return L, R
def feasible_pairs(n, pi, L, R):
    return all(R[i] & L[j] for i in range(n) for j in range(i+1, n))
# 1. literal definitions == set-pair characterisation, on random migration sets
for (k, m) in [(3,1),(3,2),(4,1),(3,3),(4,2),(5,1)]:
    n = k*m; tH, tG = times(n)
    for trial in range(400):
        pi = list(itertools.chain.from_iterable([p]*m for p in range(k))); random.shuffle(pi)
        M = sorted(candidates(n, k, pi, tH))
        S = {x for x in M if random.random() < random.choice([0.2,0.4,0.6])}
        L, R = setpairs(n, pi, S, tH)
        assert feasible_literal(n, k, pi, S, tH, tG) == feasible_pairs(n, pi, L, R), (k, m, pi, S)
print('literal coverage == set-pair condition R_i & L_j != {} on 2400 random instances')
# 2. exact MS-GC optimum by DFS over per-qubit (L_i, R_i) choices, for every balanced allocation
def opt(n, k, pi):
    opts = []
    for q in range(n):
        others = [p for p in range(k) if p != pi[q]]
        subs = [frozenset(c) for r in range(len(others)+1) for c in itertools.combinations(others, r)]
        opts.append(sorted(((frozenset({pi[q]}) | a, frozenset({pi[q]}) | b, len(a)+len(b)) for a in subs for b in subs), key=lambda x: x[2]))
    best = [10**9]
    Ls = [None]*n; Rs = [None]*n
    def dfs(q, cost):
        if cost >= best[0]: return
        if q == n: best[0] = cost; return
        for (L, R, c) in opts[q]:
            if cost + c >= best[0]: break
            if all(Rs[i] & L for i in range(q)):
                Ls[q], Rs[q] = L, R
                dfs(q+1, cost+c)
    dfs(0, 0)
    return best[0]
for (k, m) in [(3,1),(3,2),(4,1),(4,2),(5,1)]:
    n = k*m
    canon = [i//m for i in range(n)]
    oc = opt(n, k, canon)
    seen = set(); worst = None; allvals = []
    for perm in itertools.permutations(range(n)):
        pi = tuple(canon[perm[i]] for i in range(n))
        # relabel modules by first occurrence (module labels are symmetric)
        rel = {}; key = tuple(rel.setdefault(p, len(rel)) for p in pi)
        if key in seen: continue
        seen.add(key)
        o = opt(n, k, list(key)); allvals.append(o)
        assert o >= oc, (k, m, key, o, oc)
    print(f'k={k} m={m}: OPT_GC(pi*) = {oc} (MS-HC optimum m*C(k,2) = {m*k*(k-1)//2}); balanced allocations up to relabelling = {len(seen)}; min OPT_GC over them = {min(allvals)}')
print('ALL_OK')
